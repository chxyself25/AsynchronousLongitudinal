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
library(fda)
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
sd.pcs = c(8, sqrt(4), 1); K <- length(sd.pcs)
SigmaXZ <- matrix(0, ncol = dx+dz, nrow = dx+dz)
f.basis <- create.fourier.basis(rangeval = c(0, 10*dx), nbasis = ifelse(K%%2 == 0, K+1, K))
tt <- seq(0, 10, by = 0.01)
for (j in 1:dx) {
  for (l in 1:dx) {
    phij <- eval.basis(tt+(j-1)*10, f.basis)[,1:K]; phij <- sweep(phij, 2, sd.pcs, "*")
    phil <- eval.basis(tt+(l-1)*10, f.basis)[,1:K]; phil <- sweep(phil, 2, sd.pcs, "*")
    SigmaXZ[j,l] <- func_prod(phij, phil, tt)
  }
}
SigmaXZ[(dx+1):(dx+dz), (dx+1):(dx+dz)] <- diag(1, nrow = dz)

## decide fixed bandwidth in simulations 
optns <- list(dataType = "Sparse", nRegGrid = 60, methodBwMu = "GCV", methodBwCov = "GCV")
bwres <- foreach(i = 1:20, .combine = "rbind") %dopar% {
  # simulate a data set
  data.list <- sim_data2(n=n, m=m, dx=dx, dz=dz, sd.pcs = sd.pcs,
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
cat("fixed bandwithds: ", bws, "\n")
optns <- list(dataType = "Sparse", nRegGrid = 60, userBwMu = bws[1], userBwCov = bws[2])



############ results evaluation ##################
res <- readRDS("./Simulation/multix2_sim_200.rds")
xc <- res %>% 
  group_by(method) %>%
  summarize_all(function(x){mean((x))}) %>% as.data.frame
xc

