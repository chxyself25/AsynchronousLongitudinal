## multi sample splitting, assign p-value for each variable
source("../Funcs/helpers.R")
source("../Funcs/Asyncreg_funcs.R")
library(MASS)
library(Matrix, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(doParallel)
registerDoParallel(cores = 50)
library(glmnet, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
sourceDir <- function(path, trace = FALSE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../fdapace/R/")
optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "Default", methodBwCov = "Default")
dat <- readRDS(file = "./YXZ_valid_ID.rds")
IDs <- unique(dat$ID)
n <- length(IDs)
n1 <- floor(n/2)

######### resample from data #################
res <- NULL
naive <- NULL
for (b in 101:200) {
  cat("Doing ", b, "\n")
  ## first half random sample
  s.idx <- sample(IDs, n1, replace = FALSE)
  datb <- subset(dat, ID %in% s.idx)
  datb.list <- split(datb, f = as.factor(datb$ID))
  matb <- OrganizeXYZ(dat.list = datb.list, c.idx = names(dat), n = n1, ncol = ncol(datb), optns)
  MX = matb$MX
  #cat(colnames(MX), dim(MX), "\n")
  resb <- rep(1, ncol(MX)); names(resb) <- colnames(MX)
  # select variable
  lso <- cv.glmnet(x=MX, y=matb$Y, family = "gaussian", type.measure = "mse", lambda = seq(0.001, 1, length.out = 200),
            nfolds = 10, alpha = 1, intercept = FALSE, standardize = FALSE)
  beta.lso <- unlist(coef(lso, s = 'lambda.min'))
  c.idx <- colnames(MX)[which(beta.lso[-1] != 0)]
  cat("Selected: ", c.idx, "\n" )
  # use second part of data to calculate p-value
  cat("Second part of data: ", "\n")
  s2.idx <- setdiff(unique(IDs), s.idx)
  datb2 <- subset(dat, ID %in% s2.idx)
  datb2.list <- split(datb2, f = as.factor(datb2$ID))
  matb2 <- OrganizeXYZ(dat.list = datb2.list, c.idx = c.idx, n = length(datb2.list), ncol = ncol(datb2), optns)
  MX2 <- matb2$MX; Y2 <- matb2$Y
  #cat(colnames(MX2), "\n")
  beta <- solve(t(MX2)%*% MX2) %*% t(MX2) %*% Y2
  sigma2 <- (t(Y2)%*%Y2 - t(Y2) %*% MX2 %*% beta)/(length(Y2)-length(c.idx))
  test.stat <- beta/sqrt(diag(c(sigma2)*solve(t(MX2)%*%MX2)))
  prob <- 2*pt(abs(test.stat), df = length(Y2)-length(c.idx), lower.tail = FALSE) # two side p-value
  resb0 <- resb
  resb0[c.idx] <- sapply(prob*length(c.idx), function(x) {min(1, x)})
  naive <- rbind(naive, resb0)
  cat("bootstrap procedure for estimating sde: ", "\n")
  #cat(test.stat, prob, "\n")
  sde <- BSDbeta(datb2.list, c.idx, B = 250, optns)
  test.stat <- beta/sde
  prob <- 2*pnorm(abs(test.stat), mean = 0, sd = 1, lower.tail = FALSE)
  resb[c.idx] <- sapply(prob*length(c.idx), function(x) {min(1, x)})
  cat(resb, "\n")
  res <- rbind(res, resb)
  saveRDS(list(naive = naive, bsd = res), file = "./resample_fm_lasso_pvalue_sub200_b250.rds")
}
saveRDS(list(naive = naive, bsd = res), file = "./resample_fm_lasso_pvalue_sub200_b250.rds")
