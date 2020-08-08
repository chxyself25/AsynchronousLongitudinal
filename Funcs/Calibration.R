## functional match methods, two steps: first selecting variables, then estimate betas
FMlasso <- function(pca.list, t, Y, Z, n, m) {
  dx <- length(pca.list)
  Xhat <- matrix(0, ncol = dx, nrow = m*n)
  for (v in 1:dx) {
    pca <- pca.list[[v]]
    s.range <- range(pca$workGrid)
    pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + pca$xiEst %*% t(pca$phi)
    #pred <- pca$xiEst %*% t(pca$phi)
    impute_res = foreach(j = 1:n, .combine = "c") %dopar% {
      tt <- t[j,]
      in.idx <- tt >= s.range[1] & tt <= s.range[2]
      resj <- rep(NA, m)
      resj[in.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[in.idx], mu = pred[j,])
      resj
    }
    Xhat[,v] <- impute_res
  }
  # 2 step lasso
  MX <- cbind(Xhat, Z)
  #cat(ncol(MX), "\n")
  nna.idx <- apply(MX, 1, function(x) {all(!is.na(x))})
  lasso.res <- lasso_2step0(MX[nna.idx,], Y[nna.idx], nfolds = 10)
  return(list(Xhat=Xhat, betas = lasso.res$betas, selection = lasso.res$c.idx))
}

## multivariate functional match method: estimate xi conditioning on observations from all variables
MFMlasso <- function(pca_cov, s, t, Y, Z, n, m) {
  pca.list <- pca_cov$pca.list
  dx <- length(pca.list)
  ## FullPACE pc scores
  xifull <- GetFullPCScore(pca.list = pca.list, cr.cov = pca_cov$cr.cov, W_is = pca_cov$W_is, s = s, n = n, m = m)
  Xtail <- matrix(0, ncol = dx, nrow = m*n)
  for (v in 1:dx) {
    pca <- pca.list[[v]]
    pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + xifull[[v]] %*% t(pca$phi)
    #pred <- xifull[[v]] %*% t(pca$phi)
    s.range <- range(pca$workGrid)
    impute_res = foreach(j = 1:n, .combine = "c") %dopar% {
      tt <- t[j,]
      in.idx <- tt >= s.range[1] & tt <= s.range[2]
      resj <- rep(NA, m)
      resj[in.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[in.idx], mu = pred[j,])
      resj
    }
    Xtail[,v] <- impute_res
  }
  # 2 step lasso
  MX <- cbind(Xtail, Z)
  nna.idx <- apply(MX, 1, function(x) {all(!is.na(x))})
  lasso.res <- lasso_2step0(MX[nna.idx,], Y[nna.idx], nfolds = 10)
  return(list(Xtail = Xtail, betas = lasso.res$betas, selection = lasso.res$c.idx))
}

## function for last observation carried forward
LOCFlasso <- function(W_s, s, t, Y, Z, n, m) {
  dx <- ncol(W_s)
  nXhat <- matrix(0, ncol = dx, nrow = m*n)
  for (v in 1:dx) {
    Wsv <- matrix(W_s[,v], ncol = m, nrow = n, byrow = TRUE)
    sv <- s[[v]]
    nXhat[,v] <- unlist(lapply(1:n, function(i) {sapply(t[i,], function(y) {last_obs(Wsv[i,], sv[i,], y, nearest = FALSE)})}))
  }
  #nXhat <- apply(nXhat, 2, function(x) {x - mean(x, na.rm = TRUE)})
  # 2 step lasso
  MX <- cbind(nXhat, Z)
  nna.idx <- apply(MX, 1, function(x) {all(!is.na(x))})
  lasso.res <- lasso_2step0(MX[nna.idx,], Y[nna.idx], nfolds = 10)
  return(list(Xhat = nXhat, betas = lasso.res$betas, selection = lasso.res$c.idx))
}

# function that can return the nearest last observation indexes of s
# if there is no last observation, then return the nearest observation
last_obs <- function(ws=NULL, s, tt, nearest = FALSE) {
  idx1 <- which(s < tt)
  if (length(idx1) == 0) {
    if (!nearest) {return(NA)}
    else {return(ws[which.min(abs(s-tt))])}
  } else {
    return(ws[max(idx1)])
    #return(ws[idx1[which.max(s[idx1])]])
  }
}
