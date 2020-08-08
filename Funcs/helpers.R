## functions for organizing dataset (normalization, imputation), used in real data case
## functional standardize Y variable: global cognitive function score
FStandardizeYZ <- function(dat.list, n, ncol, optns = list()) {
  Lt <- lapply(dat.list, '[[', 'COGDAY')
  Ly <- lapply(dat.list, '[[', 'GLBSCORE')
  ## standardize Y variable
  nna.idx <- lapply(1:length(Ly), function(i) {intersect(which(!is.na(Lt[[i]])), which(!is.na(Ly[[i]])))})
  Lt <- lapply(1:n, function(i) {Lt[[i]][nna.idx[[i]]]})
  Ly <- lapply(1:n, function(i) {Ly[[i]][nna.idx[[i]]]})
  pca <- FPCAIC(Ly, Lt, optns, spec = list("AIC", FALSE))
  Y <- unlist(lapply(1:n, function(i) {pca$resid[[i]]/sqrt(diag(pca$fittedCovObs[[i]]) + pca$sigma2)}))
  Zmat <- do.call("rbind", lapply(1:n, function(i) {dat.list[[i]][nna.idx[[i]], 19:ncol]}))
  Zmat[,'RACE'] <- NULL
  ## standardize Z variable
  for (v in names(Zmat)) {
    Zv <- Zmat[,v]
    if (length(unique(Zv)) > 2) {
      Zmat[,v] <- (Zv - mean(Zv, na.rm = TRUE))/sd(Zv, na.rm = TRUE)
    }
  }
  return(list(Y = Y, Zmat = Zmat, Ly = Ly, Lt = Lt))
}

## imputation by FPCA for each X variable
FMimputeX <- function(dat.list, xvs, Lt, n, optns = list()) {
  XX <- NULL
  for (type in c("PHY", "HRM", "CVR")) {
    xv <- xvs[[type]]
    if (length(xv) == 0) {next} 
    Ls <- lapply(dat.list, '[[', paste0(type, "DAY"))
    ## design matrix for variable in category type
    X <- foreach (v = xv, .combine = 'cbind') %do% {
      Lw0 <- lapply(dat.list, '[[', v)
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
      impute_res
    }
    XX <- cbind(XX, X)
  }
 colnames(XX) <- unlist(xvs)
 return(XX)
}

## wrap up last two functions, standardization and imputation
OrganizeXYZ <- function(dat.list, c.idx, n, ncol, optns = list()) {
  yz.res <- FStandardizeYZ(dat.list, n, ncol, optns)
  Zmat <- yz.res$Zmat; Y <- yz.res$Y
  z.idx <- intersect(names(Zmat), c.idx)
  xvs <- list(PHY = c("BMI", "DIABP", "SYSBP", "PULSE"), HRM = c("DHAS", "SHBG", "FSH"),
            CVR = c("CHOLRES", "GLUCRES", "INSURES", "TRIGRES"))
  x.idx <- lapply(xvs, function(x) {intersect(x, c.idx)})
  XX <- FMimputeX(dat.list, xvs = x.idx, yz.res$Lt, n, optns)
  MX <- as.matrix(cbind(Zmat[,z.idx], XX))
  nna.idx <- apply(MX, 1, function(x) {all(!is.na(x))})
  MX <- MX[nna.idx,]; Y <- Y[nna.idx]
  return(list(MX = MX, Y = Y))
}

## bootstrap estimate standard error of beta
BSDbeta <- function(dat.list, c.idx, B = 250, optns = list()) {
  n <- length(dat.list)
  b.res <- foreach(b = 1:B, .combine = "rbind") %dopar% {
    b.idx <- sample(1:n, n, replace = TRUE)
    b.dat.list <- dat.list[b.idx]
    b.mat <- OrganizeXYZ(b.dat.list, c.idx, n, ncol = ncol(b.dat.list[[1]]), optns)
    b.MX <- b.mat$MX; b.Y <- b.mat$Y
    beta <- tryCatch(solve(t(b.MX)%*% b.MX) %*% t(b.MX) %*% b.Y, error = function(e) {rep(NA, ncol(b.MX))})
    c(beta)
  }
  res <- apply(b.res, 2, sd, na.rm = TRUE)
  return(res)
}



