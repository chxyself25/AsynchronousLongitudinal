## simulation for one X case
## compare with simple last observation carried forward and Cao's paper
source("./Funcs/sim_funcs.R")
source("./Funcs/Asyncreg_funcs.R")
library(fdapace)
library(MASS)
library(doParallel)
library(AsynchLong)
library(caTools)
library(doParallel)
registerDoParallel(cores = 50)
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
psi <- function(t, k) {
  if (k==1) {
    t <- as.matrix(t)
    matrix(1/sqrt(10), ncol = ncol(t), nrow = nrow(t))
  } else if (k == 2) {
    (1/sqrt(5))*sin(2*pi*t/10)
  } else {
    (1/sqrt(5))*cos(2*pi*t/10)
  }
}

# error variance structure
evar_fun <- function(tt) {
  covs <- outer(tt, tt, function(t1,t2) {
    2^(-abs(t1-t2)/5)
  })
  return(covs)
}
b <- c(1, 2); n <- 100; m <- 5; grid.num <- 60
sd.pcs <- c(2, sqrt(2), 1)

## test for simulation function
##data.list <- sim_datax(n=n, m=m, sd.pcs = sd.pcs, sd.e = 1, dx = 1, beta = b, method = "fix")

## simulation 200 times
res <- foreach(sim = 1:500, .combine = "rbind") %do% {
  ## univariate representation simulation
  data.list <- sim_datax(n=n, m=m, sd.pcs = sd.pcs, sd.e = 1, dx = 1, beta = b, method = "fix")
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
  fm.res <- FMlassoX1(Y, t, Lw, Ls, optns = NULL)
  b.fm <- fm.res$beta
  # bootstrap sampling for beta standard error estimation
  b1.bts <- foreach(bts = 1:8, .combine = "c") %dopar% {
    s.idx <- sample(1:n, size = 100, replace = TRUE)
    Lwb <- Lw[s.idx]
    Lsb <- Ls[s.idx]
    tb <- t[s.idx,]
    Yb <- c(matrix(Y, ncol = n, byrow = FALSE)[,s.idx])
    fm.b <- FMlassoX1(Yb, tb, Lwb, Lsb, optns = NULL)
    fm.b$beta[2]
  }
  # simple last observation carried forward, same result with the package
  locf.res <- LOCFlassoX(Y, t, W_s, s)
  b.locf <- locf.res$beta
  # kernel last observation carried forward
  data.x <- data.frame(ID = rep(1:n, each = m), s = c(t(s[[1]])), W = W_s)
  data.y <- data.frame(ID = rep(1:n, each = m), t = c(t(t)), Y = Y)
  ti.res <- asynchTI(data.x, data.y)
  b.ti <- ti.res$betaHat
  # summarize results
  betas <- cbind(b.fm, b.locf, as.matrix(b.ti))
  data.frame(beta0 = betas[1,], beta1 = betas[2,], method = c("fm", "locf", "ti"),
             pca.K = c(fm.res$pca.K, NA, NA), pca.sigma2 = c(fm.res$pca.sigma2, NA, NA), 
             pca.rho = c(fm.res$pca.rho, NA, NA), naive.sigma2 = c(fm.res$beta.sigma2, locf.res$beta.sigma2, NA),
             std = c(sd(b1.bts), NA, ti.res$stdErr[2]))
}
saveRDS(res, file = "./onex_sim_fix_500.rds")

########## summary simulation results ################
## settings1: increasing trend mean, n = 200, independent/dependent error variance
## settings2, n = 200, independent/dependent error variance
library(dplyr)
# true beta and covariance of X in simulation
b0 <- 1; b1 <- 2; sd.pcs <- c(2, sqrt(2), 1)
sigmaX <- sum(sd.pcs^2)
fm.res <- list()
for (set in 1:2) {
  for (ind in c("", "n")) {
    resi <- readRDS(paste0("./oneX/Simulations/onex", set, "_sim_", ind, "ind_200.rds"))
    ## evaluation of functional match method
    fm.resi <- subset(resi, method == "fm") %>%
      summarize(SE = mean((beta0-b0)^2+(beta1-b1)^2), SE.std = sd((beta0-b0)^2+(beta1-b1)^2),
                PE = mean(sigmaX*(beta1-b1)^2), PE.std = sd(sigmaX*(beta1-b1)^2),
                bias0 = mean(beta0-b0), bias0.std = sd(beta0-b0), 
                bias1 = mean(beta1-b1), bias1.std = sd(beta1-b1),
                #beta.std0 = sd(beta0), beta.std1 = sd(beta1), 
                naive0 = mean(sqrt(naive.sigma20)), naive.ci0 = mean(abs(beta0-mean(beta0)) < 1.96*sqrt(naive.sigma20)), 
                naive1 = mean(sqrt(naive.sigma21)), naive.ci1 = mean(abs(beta1-mean(beta1)) < 1.96*sqrt(naive.sigma21)),
                bts0 = mean(std0), bts.ci0 = mean(abs(beta0-mean(beta0)) < 1.96*std0), 
                bts1 = mean(std1), bts.ci1 = mean(abs(beta1-mean(beta1)) < 1.96*std1)) %>% as.data.frame
    fm.resi <- apply(fm.resi, 2, round, 4)
    # output in format
    strsi <- c()
    for (j in seq(1, 8, by = 2)) {
      strsi <- c(strsi, paste0(fm.resi[j], "(", fm.resi[j+1], ")"))
    }
    names(strsi) <- names(fm.resi)[seq(1, 8, by = 2)]
    strsi <- c(strsi, fm.resi[9:16])
    fm.res[[paste0("set", set, ind)]] <- strsi
  }
}
tbl <- as.data.frame(fm.res)
out.names <- c("SE", "PE", "$\\Delta \\beta_0$", "$\\Delta \\beta_x$", "$s.e.^*(\\wh \\beta_0)$", "$\\mathrm{CP}^*(\\wh \\beta_0)$",
               "$s.e.^*(\\wh \\beta_x)$", "$\\mathrm{CP}^*(\\wh \\beta_x)$",
               "$s.e.(\\wh \\beta_0)$", "$\\mathrm{CP}(\\wh \\beta_0)$", 
               "$s.e.(\\wh \\beta_x)$", "$\\mathrm{CP}(\\wh \\beta_x)$")
names(out.names) <- rownames(tbl)
for (i in rownames(tbl)) {
  cat(paste(out.names[i], paste(unlist(tbl[i,]), collapse = " & "), sep = " & "), "\\\\", "\n")
}

# comparing with other methods
out.names <- c("SE", "PE", "$\\Delta \\beta_0$", "$\\Delta \\beta_x$", 
               "$s.e.(\\wh \\beta_0)$", "$\\mathrm{CP}(\\wh \\beta_0)$",
               "$s.e.(\\wh \\beta_x)$", "$\\mathrm{CP}(\\wh \\beta_x)$")
for (set in 1:2) {
  cat("setting ", set, "\n")
  tbli <- list()
  for (ind in c("", "n")) {
    resi <- readRDS(paste0("./oneX/Simulations/onex", set, "_sim_", ind, "ind_200.rds"))
    all.resi <- resi %>%
      group_by(method) %>%
      summarize(SE = mean((beta0-b0)^2+(beta1-b1)^2), SE.std = sd((beta0-b0)^2+(beta1-b1)^2),
                PE = mean(sigmaX*(beta1-b1)^2), PE.std = sd(sigmaX*(beta1-b1)^2),
                bias0 = mean(beta0-b0), bias0.std = sd(beta0-b0),
                bias1 = mean(beta1-b1), bias1.std = sd(beta1-b1),
                std0 = if_else(all(method == "locf"), mean(sqrt(naive.sigma20)), mean(std0)),
                ci0 = if_else(all(method == "locf"), mean(abs(beta0-mean(beta0)) < 1.96*sqrt(naive.sigma20)), mean(abs(beta0-mean(beta0)) < 1.96*std0)),
                std1 = if_else(all(method == "locf"), mean(sqrt(naive.sigma21)), mean(std1)),
                ci1 = if_else(all(method == "locf"), mean(abs(beta1-mean(beta1)) < 1.96*sqrt(naive.sigma21)), mean(abs(beta1-mean(beta1)) < 1.96*std1))) %>% as.data.frame
    all.resi$method <- toupper(replace(as.character(all.resi$method), grep("ti", all.resi$method), "czf"))
    rownames(all.resi) <- all.resi$method; all.resi$method <- NULL
    all.resi <- round(all.resi, 3)
    for (m in c("FM", "CZF", "LOCF")) {
      strim <- c()
      for (j in seq(1, 8, by = 2)) {
        strim <- c(strim, paste0(all.resi[m, j], "(", all.resi[m, j+1], ")"))
      }
      names(strim) <- names(all.resi)[seq(1, 8, by = 2)]
      strim <- c(strim, unlist(all.resi[m, 9:12]))
      tbli[[paste0(ind, "ind", m)]] <- strim
    }
  }
  tbli <- as.data.frame(tbli)
  names(out.names) <- rownames(tbli)
  for (j in rownames(tbli)) {
    cat(paste(out.names[j], paste(unlist(tbli[j,]), collapse = " & "), sep = " & "), "\\\\", "\n")
  }
}


