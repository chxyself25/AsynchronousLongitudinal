# use impuated matrix to do lasso regression
# simulation to evaluate is this is better than last observation carry forward
library(pracma)
library(mvtnorm)
library(doParallel)
registerDoParallel(cores = 8)

#####################################################
##############Functional Match Simulation############
#####################################################
# use the functions in sim_A11.R
sizes <- c(100, 200, 300, 400, 500)
res <- list()
mean_fun <- function(tt) {
  tt + sin(tt)
}
VD <- list()
for (i in 1:10) {
  VD[[i]] <- randortho(4, type = "orthonormal")/sqrt(10)
}
# repeat 100 times, use median
betan <- foreach(i = 1:200) %dopar% {
  bx <- c(0.5,1,1.5,2,0,0,0,0,0,0)
  #bz <- c(0,0,0,3,0,0,0,0,0,0)
  bz <- c(1,0,0,0,0,0,0,0,0,0)
  C.p <- c(1,2,3,4,11)
  IC.p <- (1:20)[!(1:20) %in% C.p]
  XZ <- cbind(diag(rep(0.01, 4)), matrix(0, ncol = 6, nrow = 4))
  #XZ <- rbind(diag(rep(0.01, 10)), matrix(0, ncol = 10, 30))
  #data.list <- sim.data(n, m = 5, cor.pcs = c(0.8, 0.8, 0.8), sd.pcs = c(2, sqrt(2), 1), cov.xz = 0.1, sd.e = 1)
  #data.list <- sim.data2(n, m = 5, sd.pcs = (41-(1:40))/40, cov.XZ = XZ, sd.e = 1, basis.method = "split")
  data.list <- sim.data2(n, m = 5, sd.pcs = c(2,sqrt(3),sqrt(2),1), cov.XZ = XZ, sd.e = 1, basis.method = "random")
  #saveRDS(data.list, file = paste0("./datasets/sim_", i, ".rds"))
  #data.list <- readRDS(paste0("./datasets/sim_", i, ".rds"))
  # observed X (on S domain)
  W_s <- data.list$Ws
  # true X on t domain
  X_t <- data.list$Xt
  # time-invariant variables
  Z <- data.list$Zt
  # sigma matrix
  #tt <- seq(0, 10, by = 0.05)
  #SigmaX <- getSigmaX(xeigenfuns = cbind(psi(tt, 1), psi(tt, 2), psi(tt, 3)), xeigenv = c(4,2,1),
  #xxeigenv = data.list$V[1:30, 1:30], xzeigenv = data.list$V[1:30, 31:40], dx = 10, dz = 10, t = tt)
  # observed y response
  Y <- cbind(X_t, Z) %*% c(bx, bz) + rnorm(n*5, mean = 0, sd = 1)
  # two observation domain
  s <- data.list$s
  t <- data.list$t
  ## funcional match
  fm <- FMlasso(W_s, s, t, Z=Z, Y, type = "mse", nfolds = 10)
  beta.fm <- fm$betas
  var.fm <- fm$selection
  #lasso.res <- lasso_2step(cbind(X_t, Z), Y, type = "mse", pars.n = 20, nfolds = 10, intercept = TRUE)
  #beta.fm <- lasso.res$betas
  #var.fm <- lasso.res$c.idx
  ## Multivariate functional match
  #mfm <- MFMlasso(W_s, s, t, Z=Z, Y, type = "mse", nfolds = 10, pca.list = fm$pca.list)
  #beta.mfm <- mfm$betas
  #var.mfm <- mfm$selection
  # last observation carry forward
  locf <- LOCFlasso(W_s, s, t, Z, Y, type = "mse", nfolds = 10)
  beta.locf <- locf$betas
  var.locf <- locf$selection
  # C: number of correct predictor identified
  C.num <- c(sum(C.p %in% var.fm), sum(C.p %in% var.locf)) 
  # IC: number of incorrect predictor identified
  IC.num <- c(sum(IC.p %in% var.fm), sum(IC.p %in% var.locf))
  # SE of beta
  beta <- c(0, bx, bz)
  SE <- c(sum((beta.fm-beta)^2), sum((beta.locf-beta)^2))
  # PE 
  #PE <- c(t(beta-beta.mfm)[-1] %*% SigmaX %*% (beta-beta.mfm)[-1], t(beta-beta.fm)[-1] %*% SigmaX %*% (beta-beta.fm)[-1],
  #t(beta-beta.locf)[-1] %*% SigmaX %*% (beta-beta.locf)[-1])
  bx1 <- c(beta.fm[2]-0.5, beta.locf[2]-0.5)
  bx2 <- c(beta.fm[3]-1, beta.locf[3]-1)
  bx3 <- c(beta.fm[4]-1.5, beta.locf[4]-1.5)
  bx4 <- c(beta.fm[5]-2, beta.locf[5]-2)
  bz <- c(beta.fm[12]-1, beta.locf[12]-1)
  # correct predictor estimate bias
  # data.frame(f.Cnum=C.num[1], l.Cnum=C.num[2], f.ICnum=IC.num[1], l.ICnum=IC.num[2],
  #            f.SE=SE[1], l.SE=SE[2], f.PE=PE[1], l.PE = PE[2],
  #            f.bx1= beta.fm[1]-0.5, f.bx2 = beta.fm[2]-1, f.bx3=beta.fm[3]-1.5, f.bz1=beta.fm[11]-1,
  #            l.bx1=beta.locf[1]-0.5, l.bx2 = beta.locf[2]-1, l.bx3 = beta.locf[3]-1.5, l.bz1=beta.locf[11]-1)
  #data.frame(C.num = C.num, IC.num = IC.num, SE = SE, PE = PE, bx1 = bx1, bx2 = bx2, bx3 = bx3, bx4 = bx4, bz = bz)
  data.frame(C.num = C.num, IC.num = IC.num, SE = SE, bx1 = bx1, bx2 = bx2, bx3 = bx3, bx4 = bx4, bz = bz)
}

mfm.res <- rbindlist(lapply(betan, function(x) {x[1,]}))
colMeans(mfm.res)
fm.res <- rbindlist(lapply(betan, function(x) {x[2,]}))
colMeans(fm.res)
locf.res <- rbindlist(lapply(betan, function(x) {x[3,]}))
colMeans(locf.res)
