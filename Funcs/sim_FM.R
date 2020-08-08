# use impuated matrix to do lasso regression
# simulation to evaluate is this is better than last observation carry forward
library(fdapace)
library(glmnet)
library(dplyr)
library(ggplot2)
library(MASS)
library(doParallel)
registerDoParallel(cores = 16)

#####################################################
##############Functional Match Simulation############
#####################################################
# use the functions in sim_A11.R
sizes <- c(100, 200, 300, 400, 500)
res <- list()
mean_fun <- function(tt) {
  tt + sin(tt)
}
# repeat 100 times, use median
betan <- foreach(i = 1:200) %dopar% {
  bx <- c(0.5,1,1.5,2,0,0,0,0,0,0)
  #bz <- c(0,0,0,3,0,0,0,0,0,0)
  bz <- c(1,0,0,0,0,0,0,0,0,0)
  C.p <- c(1,2,3,4,11)
  IC.p <- (1:20)[!(1:20) %in% C.p]
  data.list <- sim.data(n, m = 5, cor.pcs = c(0.6, 0.4, 0.2), sd.pcs = c(2, sqrt(2), 1), cov.xz = 0.1, sd.e = 1.5)
  # observed X (on S domain)
  W_s <- data.list$Ws
  # true X on t domain
  X_t <- data.list$Xt
  # time-invariant variables
  Z <- data.list$Zt
  # sigma matrix
  tt <- seq(0, 10, by = 0.05)
  SigmaX <- getSigmaX(xeigenfuns = cbind(psi(tt, 1), psi(tt, 2), psi(tt, 3)), xeigenv = c(4,2,1),
                      xxeigenv = data.list$V[1:30, 1:30], xzeigenv = data.list$V[1:30, 31:40], dx = 10, dz = 10, t = tt)
  # observed y response
  Y <- cbind(X_t, Z) %*% c(bx, bz) + rnorm(n*5, mean = 0, sd = 1)
  # two observation domain
  s <- data.list$s
  t <- data.list$t
  ## funcional match
  fm <- FMlasso(W_s, s, t, Z=Z, Y, type = "mse", nfolds = 10)
  beta.fm <- fm$betas
  var.fm <- fm$selection
  ## Multivariate functional match
  mfm <- MFMlasso(W_s, s, t, Z=Z, Y, type = "mse", nfolds = 10, pca.list = fm$pca.list)
  beta.mfm <- mfm$betas
  var.mfm <- mfm$selection
  # last observation carry forward
  locf <- LOCFlasso(W_s, s, t, Z, Y, type = "mse", nfolds = 10)
  beta.locf <- locf$betas
  var.locf <- locf$selection
  # C: number of correct predictor identified
  C.num <- c(sum(C.p %in% var.mfm), sum(C.p %in% var.fm), sum(C.p %in% var.locf)) 
  # IC: number of incorrect predictor identified
  IC.num <- c(sum(IC.p %in% var.mfm), sum(IC.p %in% var.fm), sum(IC.p %in% var.locf))
  # SE of beta
  beta <- c(0, bx, bz)
  SE <- c(sum((beta.mfm-beta)^2), sum((beta.fm-beta)^2), sum((beta.locf-beta)^2))
  # PE 
  PE <- c(t(beta-beta.mfm)[-1] %*% SigmaX %*% (beta-beta.mfm)[-1], t(beta-beta.fm)[-1] %*% SigmaX %*% (beta-beta.fm)[-1], 
          t(beta-beta.locf)[-1] %*% SigmaX %*% (beta-beta.locf)[-1])
  bx1 <- c(beta.mfm[2]-2.5, beta.fm[2]-0.5, beta.locf[2]-0.5)
  bx2 <- c(beta.mfm[3]-1, beta.fm[3]-1, beta.locf[3]-1)
  bx3 <- c(beta.mfm[4]-1.5, beta.fm[4]-1.5, beta.locf[4]-1.5)
  bx4 <- c(beta.mfm[5]-2, beta.fm[5]-2, beta.locf[5]-2)
  bz <- c(beta.mfm[12]-1, beta.fm[12]-1, beta.locf[12]-1)
  # correct predictor estimate bias
  # data.frame(f.Cnum=C.num[1], l.Cnum=C.num[2], f.ICnum=IC.num[1], l.ICnum=IC.num[2],
  #            f.SE=SE[1], l.SE=SE[2], f.PE=PE[1], l.PE = PE[2], 
  #            f.bx1= beta.fm[1]-0.5, f.bx2 = beta.fm[2]-1, f.bx3=beta.fm[3]-1.5, f.bz1=beta.fm[11]-1,
  #            l.bx1=beta.locf[1]-0.5, l.bx2 = beta.locf[2]-1, l.bx3 = beta.locf[3]-1.5, l.bz1=beta.locf[11]-1)
  data.frame(C.num = C.num, IC.num = IC.num, SE = SE, PE = PE, bx1 = bx1, bx2 = bx2, bx3 = bx3, bx4 = bx4, bz = bz)
}

mfm.res <- rbindlist(lapply(betan, function(x) {x[1,]}))
colMeans(mfm.res)
fm.res <- rbindlist(lapply(betan, function(x) {x[2,]}))
colMeans(fm.res)
locf.res <- rbindlist(lapply(betan, function(x) {x[3,]}))
colMeans(locf.res)

xc1 <- res[[100]] # tau = 0.75
xc2 <- res[[100]] # tau = 1
xc3 <- res[[100]] # tau = 1.25
# function for getting standard error estimation by bootstrap
getSimSD <- function(res, method = NULL) {
  sdn <- foreach(i = 1:ncol(res), .combine = 'c') %do% {
    bsi <- foreach(j = 1:500, .combine = 'c') %dopar% {
      if (method == 'median') {
        median(sample(res[,i], size = nrow(res), replace = TRUE))
      } else {
        mean(sample(res[,i], size = nrow(res), replace = TRUE))
      }
    }
    sd(bsi)
  }
  names(sdn) <- names(res)
  return(sdn)
}
res75 <- readRDS("sim_FMvsLOCF_75.rds")
res100 <- readRDS("sim_FMvsLOCF_100.rds")
res125 <- readRDS("sim_FMvsLOCF_125.rds")

round(apply(res75, 2, median), 2)
round(getSimSD(res75, method = 'median'), 2)
round(apply(res100, 2, median), 2)
round(getSimSD(res100, method = 'median'), 2)
round(apply(res125, 2, median), 2)
round(getSimSD(res125, method = 'median'), 2)
round(apply(res75, 2, sd), 2)
round(apply(res100, 2, sd), 2)
round(apply(res125, 2, sd), 2)

round(apply(res75, 2, mean), 2)
round(getSimSD(res75, method = 'mean'), 2)
round(apply(res100, 2, mean), 2)
round(getSimSD(res100, method = 'mean'), 2)
round(apply(res125, 2, mean), 2)
round(getSimSD(res125, method = 'mean'), 2)

