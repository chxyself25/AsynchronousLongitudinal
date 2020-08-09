# subsample from original data and then do imputation and analysis on it
source("../Funcs/sim_funcs.R")
source("../Funcs/Asyncreg_funcs.R")
library(MASS)
library(Matrix, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
#library(glmnet, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
#library(Matrix, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
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
optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "GCV", methodBwCov = "GCV")
n <- 350

######### resample from data #################
dat <- readRDS(file = "./YXZ_valid_ID.rds")
#dat.list <- split(dat, f = as.factor(dat$ID))
IDs <- unique(dat$ID)
#mubw <- covbw <- NA
res <- foreach(b = 1:300) %dopar% {
  s.idx <- sample(IDs, n, replace = FALSE)
  datb <- subset(dat, ID %in% s.idx)
  datb.list <- split(datb, f = as.factor(datb$ID))
  Lt <- lapply(datb.list, '[[', 'COGDAY')
  Ly <- lapply(datb.list, '[[', 'GLBSCORE')
  nna.idx <- lapply(1:length(Ly), function(i) {intersect(which(!is.na(Lt[[i]])), which(!is.na(Ly[[i]])))})
  Lt <- lapply(1:n, function(i) {Lt[[i]][nna.idx[[i]]]})
  Ly <- lapply(1:n, function(i) {Ly[[i]][nna.idx[[i]]]})
  cat("standardization for Y", "\n")
  pca <- FPCAIC(Ly, Lt, optns, spec = list("AIC", FALSE))
  bw.res <- list("Y" = c(pca$bwMu, pca$bwCov))
  Y <- unlist(lapply(1:n, function(i) {pca$resid[[i]]/sqrt(diag(pca$fittedCovObs[[i]]) + pca$sigma2)}))
  Zmat <- do.call("rbind", lapply(1:n, function(i) {datb.list[[i]][nna.idx[[i]], 19:ncol(datb)]}))
  Zmat[,'RACE'] <- NULL
  ## standardize Z variable
  for (v in names(Zmat)) {
  Zv <- Zmat[,v]
  if (length(unique(Zv)) > 2) {
    Zmat[,v] <- (Zv - mean(Zv, na.rm = TRUE))/sd(Zv, na.rm = TRUE)
   }
  }
  xvs <- list(PHY = c("BMI", "DIABP", "SYSBP", "PULSE"), HRM = c("DHAS", "SHBG", "FSH"),
            CVR = c("CHOLRES", "GLUCRES", "INSURES", "TRIGRES"))
  # imputation for each individual variable
  cat("start doing imputation for X", "\n")
  XX <- Zmat
  for (type in c("PHY", "HRM", "CVR")) {
   #cat(type, "\n")
    xv <- xvs[[type]]
    Ls <- lapply(datb.list, '[[', paste0(type, "DAY"))
    #cat(range(unlist(Ls), na.rm = TRUE), "\n")
    ## design matrix for physical measures
    X <- foreach (v = xv, .combine = 'cbind') %do% {
      Lw0 <- lapply(datb.list, '[[', v)
      pca <- FPCAIC(Lw0, Ls, optns, spec = list("AIC", FALSE))
      use.Ls <- pca$inputData$Lt
      Lw <- lapply(1:n, function(i) {
        resi <- rep(NA, length(Ls[[i]]))
        resi[match(use.Ls[[i]], Ls[[i]])] <- pca$resid[[i]]/sqrt(diag(pca$fittedCovObs[[i]]) + pca$sigma2)
      })
      pred <- sweep(pca$xiEst %*% t(pca$phi), 2, sqrt(diag(pca$fittedCov)+pca$sigma2), "/")
      #pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = 51) + pca$xiEst %*% t(pca$phi)
      s.range <- range(pca$workGrid)
      impute_res = foreach(j = 1:n, .combine = "c") %do% {
        tt <- Lt[[j]]; ss <- Ls[[j]]
        # see if any already matched
        resj <- sapply(tt, function(x) {if (x %in% ss) {Lw[[j]][which(ss == x)]} else {NA}})
        in.idx <- which(tt >= s.range[1] & tt <= s.range[2])
        na.idx <- intersect(which(is.na(resj)), in.idx)
        if (length(na.idx) != 0) {
          resj[na.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[na.idx], mu = pred[j,])
        }
        resj
      }
      c(c(pca$bwMu, pca$bwCov), impute_res)
    }
    XX <- cbind(XX, X[3:nrow(X),])
    for (j in 1:length(xv)) {
     bw.res[[xv[j]]] <- X[1:2, j]
   }
  }
  colnames(XX)[(which(colnames(XX) == 'INSUR')+1):ncol(XX)] <- unlist(xvs)
  XX <- as.matrix(XX)
  # do lasso using imputed matrix
  cat("lasso", "\n")
  nna.idx <- apply(XX, 1, function(x) {all(!is.na(x))})
  XX <- XX[nna.idx,]; Y <- Y[nna.idx]
  lassob <- glmnet(XX, Y, family = "gaussian", lambda = seq(0.001,0.3,length.out = 200), alpha = 1, intercept = FALSE, standardize = FALSE)
  list(bw.res, rbind(lassob$lambda, lassob$beta))
}
res1 <- lapply(res, '[[', 1)
saveRDS(res1, file = "./resample_fm_bw_b300_n350.rds")
res2 <- lapply(res, '[[', 2)
saveRDS(res2, file= "./resample_fm_lasso_b300_n350_GCV.rds")
