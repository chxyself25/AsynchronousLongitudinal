## simulation for one X case
## compare with simple last observation carried forward and Cao's paper
source("./Funcs/sim_funcs.R")
source("./Funcs/Asyncreg_X1.R")
library(fdapace)
library(MASS)
library(doParallel)
library(AsynchLong)
library(caTools)
library(tvReg)
library(doParallel)
registerDoParallel(cores = 8)
sourceDir <- function(path, trace = FALSE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("./fdapace/R/")

## specification
# mean function
mean_fun <- function(tt) {
  tt + sin(tt)
  #sin(tt)
}
# define eigen functions: k = 1,2,3
psi <- function(t, k) {
  (1/sqrt(5))*sin(pi*k*t/10)
}
# error variance structure
evar_fun <- function(tt, ind = FALSE) {
  if (ind) {
    covs <- diag(1.5, nrow = length(tt))
  } else {
    covs <- outer(tt, tt, function(t1,t2) {
      2^(-abs(t1-t2)/5)
    })}
  return(covs)
}

beta_fun <- function(tt, intercept = FALSE) {
  if (intercept) {
    betat <- 0.4*tt + 0.5
  } else {
    betat <- sin(2*pi*tt/10)
  }
  return(betat)
}

n <- 200; m <- 5
sd.pcs <- c(2, sqrt(2), 1)
optns <- list(dataType = "Sparse", nRegGrid = 60, methodBwMu = "Default", methodBwCov = "Default")

## simulation 200 times
res <- foreach(sim = 1:100, .combine = "rbind") %do% {
  ## univariate representation simulation
  data.list <- sim_dataxt1(n=n, m=m, uni.f = 1, sd.pcs = sd.pcs, sd.e = 1)
  # observed X (on S domain)
  W_s <- data.list$Ws
  # true X on t domain
  X_t <- data.list$Xt
  # two observed sampling time points
  s <- data.list$s
  t <- data.list$t
  # observed y response
  Y <- data.list$Y
  
  # functional match 
  Lw <- split(W_s, as.factor(rep(1:n, each = m)))
  Ls <- lapply(1:n, function(x) {s[[1]][x,]})
  #Ly <- split(Y, as.factor(rep(1:n, each = m)))
  fm.res <- tvFMlassoX1(Y, t, Lw, Ls, optns)
  # simple last observation carried forward, same result with the package
  locf.res <- tvLOCFlassoX1(Y, t, Lw, Ls)
  # kernel last observation carried forward
  data.x <- data.frame(ID = rep(1:n, each = m), s = c(t(s[[1]])), W = W_s)
  data.y <- data.frame(ID = rep(1:n, each = m), t = c(t(t)), Y = Y)
  timepts <- fm.res$time
  #timepts <- seq(0.1, 10, by = 0.9)
  ti.res <- foreach(t = timepts, .combine = "rbind") %do% {
    timodel <- tryCatch(asynchTD(data.x, data.y, times = t, kType = "epan", lType = "identity",
             bw = fm.res$bw, verbose = FALSE), error = function(e) {NA})
    if (is.na(timodel)) {
      rep(NA, 2)
    } else {
      timodel$betaHat
    }
  }
  diff.mat <- cbind(beta_fun(timepts, intercept = TRUE), beta_fun(timepts, intercept = FALSE)) - ti.res
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {
    nna.idx <- which(!is.na(x)); trapzRcpp(timepts[nna.idx], abs(x[nna.idx]))})
  mseint <- apply(diff.mat, 2, function(x) {
    nna.idx <- which(!is.na(x)); trapzRcpp(timepts[nna.idx], (x[nna.idx])^2)})
  #b.ti <- ti.res$betaHat
  # summarize results
  resdf <- NULL
  for (i in c("mae", "mse", "maeint", "mseint")) {
    dfm <- data.frame(value = c(fm.res[[i]], locf.res[[i]]), metric = rep(paste(i, c(0,1), sep = ""),2), 
                      method = rep(c("fm", "locf"), each = 2))
    resdf <- rbind(resdf, dfm)
  }
  resdf
}
saveRDS(res, file = "./OneX/Simulations/onext1_sim_nind_ll_200.rds")
######### summarize ################
res1 <- readRDS("./OneX/Simulations/onext1_sim_nind_ll_100.rds")
res2 <- readRDS("./OneX/Simulations/onext1_sim_nind_ll_200.rds")
res3 <- readRDS("./OneX/Simulations/onext1_sim_nind_czf_200.rds")
xc <- subset(rbind(res1, res2, res3), !grepl("int", metric))
xc$method <- gsub("FM", "fImpute", toupper(xc$method))
xc$metric <- toupper(xc$metric)
ggplot(xc, aes(method, value)) + geom_boxplot() + facet_wrap(.~metric, nrow = 2, ncol = 2, scale = "free_y") + 
  theme_bw() + ylim(c(0,20))
ggsave(paste0("./OneX/fm_locf_betat_sim_ll_boxplot.pdf"), width = 6, height = 5)
xc <- subset(rbind(res1, res2, res3), grepl("int", metric))
xc$method <- gsub("FM", "fImpute", toupper(xc$method))
xc$metric <- gsub("MAEINT", "IMAE", toupper(xc$metric))
xc$metric <- gsub("MSEINT", "IMSE", xc$metric)
ggplot(xc, aes(method, value)) + geom_boxplot() + facet_wrap(.~metric, nrow = 2, ncol = 2, scale = "free_y") + theme_bw() 
ggsave(paste0("./OneX/fm_locf_betat_sim_ll_int_boxplot.pdf"), width = 6, height = 5)


########### evaluation ####################
str <- "1m7"
res <- readRDS(paste0("./OneX/Simulations/onext1_sim_nind_", str, "_200.rds")) # m = 7, and same observations for all i
res$true_beta0 <- beta_fun(res$time, intercept = TRUE)
res$true_beta1 <- beta_fun(res$time, intercept = FALSE)
## MADE and WASE
df <- res %>% group_by(method, sim) %>% 
  summarize(MAE0 = mean(abs(beta0 - true_beta0), na.rm = TRUE),
            MAE1 = mean(abs(beta1 - true_beta1), na.rm = TRUE),
            MSE0 = mean((beta0 - true_beta0)^2, na.rm = TRUE),
            MSE1 = mean((beta1 - true_beta1)^2, na.rm = TRUE)) %>% as.data.frame
df <- melt(df, id = c("method", "sim"), variable.name = "metric", value.name = "value")
ggplot(df, aes(method, value)) + geom_boxplot() + facet_wrap(.~metric, nrow = 2, ncol = 2, scale = "free_y") + theme_bw() 
ggsave(paste0("./OneX/fm_locf_betat_sim_", str, "_boxplot.pdf"), width = 6, height = 5)

## visualize coefficient function
library(gridExtra)
res1 <- subset(res, method == "fm")
b0 <- ggplot(res1) + geom_line(aes(time, beta0), color = "black") + 
  geom_line(aes(time, true_beta0), color = "red") + theme_bw()
b1 <- ggplot(res1) + geom_line(aes(time, beta1), color = "black") + 
  geom_line(aes(time, true_beta1), color = "red") + theme_bw()
ggsave(grid.arrange(b0, b1), file = paste0("./OneX/fm_betat_sim_", str, "_plot.pdf"), width = 8, height = 5)
res2 <- subset(res, method == "locf")
b0 <- ggplot(res2) + geom_line(aes(time, beta0), color = "black") + 
  geom_line(aes(time, true_beta0), color = "red") + theme_bw()
b1 <- ggplot(res2) + geom_line(aes(time, beta1), color = "black") + 
  geom_line(aes(time, true_beta1), color = "red") + theme_bw()
ggsave(grid.arrange(b0, b1), file = paste0("./OneX/locf_betat_sim_", str, "_plot.pdf"), width = 8, height = 5)



