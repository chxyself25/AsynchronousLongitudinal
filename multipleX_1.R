# simulation of multiple linear regression
## simulation using FullPACE as the estimate
source("../Funcs/sim_funcs.R")
source("../Funcs/Calibration.R")
source("../Funcs/Lasso.R")
source("./Simulation/GetFullPACE.R")
#library(pracma)
#library(mvtnorm)
library(MASS)
library(Matrix)
library(Rcpp)
library(fdapace)
sourceDir <- function(path, trace = FALSE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../fdapace/R/")
library(glmnet)
library(doParallel)
registerDoParallel(cores = 8)

## specification
# mean function
mean_fun <- function(tt) {
  tt+sin(tt)
  #matrix(0, nrow = nrow(tt), ncol = ncol(tt))
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
    covs <- outer(tt, tt, function(t1, t2) {
      2^(-abs(t1-t2)/5)
    })}
  return(covs)
}

n <- 200; m <- 5
bx <- c(1,2,0,0,0,0,0,0,0,0); dx = 10
bz <- c(0.5,1.5,0,0,0,0,0,0,0,0); dz = 10
C.p <- c(1,2,11,12)
IC.p <- (1:20)[!(1:20) %in% C.p]
cor.pcs = c(0.8, 0.6, 0.4); sd.pcs = c(2, sqrt(2), 1)
## setting 1 SigmaX
SigmaXZ <- matrix(0, ncol = dx+dz, nrow = dx+dz)
for (j in 1:dx) {
  for (l in 1:dx) {
    SigmaXZ[j,l] <- sum((sd.pcs^2)*(cor.pcs^(abs(j-l))))
  }
} 
SigmaXZ[(dx+1):(dx+dz), (dx+1):(dx+dz)] <- diag(1, nrow = dz)

## decide fixed bandwidth in simulations 
optns <- list(dataType = "Sparse", nRegGrid = 60, methodBwMu = "GCV", methodBwCov = "GCV")
bwres <- foreach(i = 1:20, .combine = "rbind") %dopar% {
  # simulate a data set
  data.list <- sim_data1(n=n, m=m, dx=dx, dz=dz, cor.pcs=cor.pcs, sd.pcs = sd.pcs, 
                         sd.e = 1, beta = c(bx,bz))
  # observed X (on S domain)
  Ws <- data.list$Ws
  # two observation domain
  s <- data.list$s
  bwsi = foreach(v = 1:dx, .combine = "rbind") %dopar% {
    # do FPCA, impute W_s on T domain 
    Ly <- split(Ws[, v], as.factor(rep(1:n, each = m)))
    Lt <- lapply(1:n, function(x) {s[[v]][x,]})
    pca <- FPCAIC(Ly, Lt, optns, spec = list("AIC", FALSE))
    c(pca$bwMu, pca$bwCov)
  }
  bwsi
}
bws <- colMeans(bwres)
optns <- list(dataType = "Sparse", nRegGrid = 60, userBwMu = bws[1], userBwCov = bws[2])
## simulate 200 times 
res <- foreach(i = 1:200, .combine = "rbind") %dopar% {
  # simulate a data set
  data.list <- sim_data1(n=n, m=m, dx=dx, dz=dz, cor.pcs=cor.pcs, sd.pcs = sd.pcs, 
                         sd.e = 1, beta = c(bx,bz))
  # observed X (on S domain)
  Ws <- data.list$Ws
  # time invariant variable Z
  Z <- data.list$Z
  # two observation domain
  s <- data.list$s
  t <- data.list$t
  # observed y response
  Y <- data.list$Y
  ypca <- FPCAIC(split(Y, as.factor(rep(1:n, each = m))), lapply(1:n, function(x) {t[x,]}), optns, spec = list("AIC", FALSE))
  Y <- unlist(ypca$resid)
  
  ## univariate FPCA and covariance calculation
  pca_cov <- uFPCACov(W_s = Ws, s = s, n = n, m = m, optns)
  pca.list <- pca_cov$pca.list
  
  ## Univariate functional match
  lassoh <- FMlasso(pca.list, t, Y, Z, n, m)
  dfh <- data.frame(C.num = sum(C.p %in% lassoh$selection), IC.num = sum(IC.p %in% lassoh$selection),
                    SE = sum((lassoh$betas - c(bx,bz))^2), PE = t(lassoh$betas-c(bx,bz))%*%SigmaXZ%*%(lassoh$betas-c(bx,bz)),
                    bx1 = lassoh$betas[1]-bx[1], bx2 = lassoh$betas[2]-bx[2],
                    bz1 = lassoh$betas[11]-bz[1], bz2 = lassoh$betas[12]-bz[2])
  
  ## Multivariate Functional Match
  lassot <- MFMlasso(pca_cov, s, t, Y, Z, n, m)
  dft <- data.frame(C.num = sum(C.p %in% lassot$selection), IC.num = sum(IC.p %in% lassot$selection),
                    SE = sum((lassot$betas - c(bx,bz))^2), PE = t(lassot$betas-c(bx,bz))%*%SigmaXZ%*%(lassot$betas-c(bx,bz)),
                    bx1 = lassot$betas[1]-bx[1], bx2 = lassot$betas[2]-bx[2],
                    bz1 = lassot$betas[11]-bz[1], bz2 = lassot$betas[12]-bz[2])
  ## last observation carried forward
  lassoc <- LOCFlasso(Ws, s, t, Y, Z, n, m)
  dfc <- data.frame(C.num = sum(C.p %in% lassoc$selection), IC.num = sum(IC.p %in% lassoc$selection),
                    SE = sum((lassoc$betas - c(bx,bz))^2), PE = t(lassoc$betas-c(bx,bz))%*%SigmaXZ%*%(lassoc$betas-c(bx,bz)),
                    bx1 = lassoc$betas[1]-bx[1], bx2 = lassoc$betas[2]-bx[2],
                    bz1 = lassoc$betas[11]-bz[1], bz2 = lassoc$betas[12]-bz[2])
  resi <- rbind(dft, dfh, dfc)
  resi$method <- c("FM", "MFM", "LOCF")
  resi
}

############ results evaluation ##################
res <- readRDS("./MultipleX/Simulation/multix2_t+sint_sim_200.rds")
xc <- res %>% 
  group_by(method) %>%
  summarize_all(function(x){round(mean(x), 4)}) %>% as.data.frame
paste(apply(xc, 1, function(x) {paste(x, collapse = " & ")}), collapse = " \\ ")

