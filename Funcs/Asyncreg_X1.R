## asynchronous longitudinal regression for one X case
## wrap up FPCA, calibration, and regression processses
## input: two functional observations including both Y and t in a list format
FMlassoX1 <- function(Y, t, Lw, Ls, optns = NULL) {
  if (is.null(optns)) {
    optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "GCV", methodBwCov = "GCV")
  }
  n <- length(Lw) # Y and W should have the same number of subjects
  pca <- FPCAIC(Lw, Ls, optns, spec = list("AIC", FALSE))
  pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + pca$xiEst %*% t(pca$phi)
  s.range <- range(pca$workGrid)
  impute_res = foreach(j = 1:n, .combine = "c") %dopar% {
    tt <- t[j,]
    in.idx <- tt >= s.range[1] & tt <= s.range[2]
    resj <- rep(NA, length(tt))
    resj[in.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[in.idx], mu = pred[j,])
    resj
  }
  nna.idx <- which(!is.na(impute_res))
  X.fm <- cbind(1, impute_res[nna.idx])
  Y.nna <- Y[nna.idx]
  b.fm <- solve(t(X.fm) %*% X.fm) %*% t(X.fm) %*% Y.nna
  # test on beta_x
  sigma2 <- (t(Y.nna)%*%Y.nna - t(Y.nna)%*% X.fm %*% b.fm)/(length(Y.nna)-2)
  #cat(t(X.fm)%*%X.fm, "\n")
  beta.sigma2 <- c(sigma2) *  solve(t(X.fm)%*%X.fm)
  test.stat <- b.fm[2]/sqrt(beta.sigma2[2,2])
  prob <- 2*pt(abs(test.stat), df = length(Y.nna)-2, lower.tail = FALSE)
  return(list(beta = b.fm, beta.sigma2 = diag(beta.sigma2),tval = test.stat, pval = prob,
              pca.sigma2 = pca$sigma2, pca.rho = pca$rho, pca.K = pca$selectK, pca.mubw = pca$bwMu, pca.covbw = pca$bwCov))
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

## function for last observation carried forward
LOCFlassoX <- function(Y, t, Ws, s) {
  n <- nrow(t)
  m <- ncol(t)
  dx <- ncol(Ws)
  nXhat <- matrix(0, ncol = dx, nrow = m*n)
  for (v in 1:dx) {
    Wsv <- matrix(Ws[,v], ncol = m, nrow = n, byrow = TRUE)
    sv <- s[[v]]
    nXhat[,v] <- c(sapply(1:n, function(i) {
      sapply(t[i,], function(y) {last_obs(Wsv[i,], sv[i,], y, nearest = FALSE)})
      }))
  }
  nna.idx <- which(apply(nXhat, 1, function(x) {all(!is.na(x))}))
  Xmat <- cbind(1, nXhat[nna.idx,])
  Y.nna <- Y[nna.idx]
  beta <- solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% Y.nna
  # test on beta_x
  sigma2 <- (t(Y.nna)%*%Y.nna - t(Y.nna) %*% Xmat %*% beta)/(length(Y.nna)-2)
  beta.sigma2 <- c(sigma2) * solve(t(Xmat)%*%Xmat)
  test.stat <- beta[2]/sqrt(beta.sigma2[2,2])
  prob <- 2*pt(abs(test.stat), df = length(Y.nna)-2, lower.tail = FALSE)
  return(list(beta = beta, beta.sigma2 = diag(beta.sigma2), tval = test.stat, pval = prob))
}

# time-dependent coefficient: only one X
tvFMlassoX1 <- function(Y, t, Lw, Ls, optns = NULL) {
  if (is.null(optns)) {
    optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "GCV", methodBwCov = "GCV")
  }
  n <- length(Lw) # Y and W should have the same number of subjects
  pca <- FPCAIC(Lw, Ls, optns, spec = list("AIC", FALSE))
  pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + pca$xiEst %*% t(pca$phi)
  s.range <- range(pca$workGrid)
  impute_res = foreach(j = 1:n, .combine = "c") %dopar% {
    tt <- t[j,]
    in.idx <- tt >= s.range[1] & tt <= s.range[2]
    resj <- rep(NA, length(tt))
    resj[in.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[in.idx], mu = pred[j,])
    resj
  }
  nna.idx <- which(!is.na(impute_res))
  tt <- c(t(t))[nna.idx]
  df <- data.frame(y = Y[nna.idx][order(tt)], x = impute_res[nna.idx][order(tt)])
  tvmodel <- tvLM(y~x, data = df, cv.block = n/10, est="ll")
  ttsort <- tt[order(tt)]
  diff.mat <- cbind(beta_fun(ttsort, intercept = TRUE), beta_fun(ttsort, intercept = FALSE)) - tvmodel$coefficients
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, abs(x))})
  mseint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, x^2)})
  return(list(mae = mae, mse = mse, maeint = maeint, mseint = mseint, coef= tvmodel$coefficients, time = ttsort,
              bw = tvmodel$bw, pca.K = pca$selectK, pca.mubw = pca$bwMu, pca.covbw = pca$bwCov))
}

tvLOCFlassoX1 <- function(Y, t, Lw, Ls) {
  n <- nrow(t)
  m <- ncol(t)
  impute_res <- foreach(j = 1:n, .combine = "c") %dopar% {
    tt <- t[j,]
    sapply(tt, function(x) {last_obs(Lw[[j]], Ls[[j]], x, nearest = FALSE)})
  }
  nna.idx <- which(!is.na(impute_res))
  tt <- c(t(t))[nna.idx]
  df <- data.frame(y = Y[nna.idx][order(tt)], x = impute_res[nna.idx][order(tt)])
  tvmodel <- tvLM(y~x, data = df, cv.block = n/10, est="ll")
  ttsort <- tt[order(tt)]
  diff.mat <- cbind(beta_fun(ttsort, intercept = TRUE), beta_fun(ttsort, intercept = FALSE)) - tvmodel$coefficients
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, abs(x))})
  mseint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, x^2)})
  return(list(mae = mae, mse = mse, maeint = maeint, mseint = mseint, 
              bw = tvmodel$bw, coef = tvmodel$coefficients, time = ttsort))
}
############# old functions using Fan's 2 step method ###############
## time-dependent coefficient: only one X 
# tdFMlassoX1 <- function(Ly, t, Lw, Ls, optns = NULL) {
#   if (is.null(optns)) {
#     optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "GCV", methodBwCov = "GCV")
#   }
#   n <- length(Lw) # Y and W should have the same number of subjects
#   pca <- FPCAIC(Lw, Ls, optns, spec = list("AIC", FALSE))
#   pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + pca$xiEst %*% t(pca$phi)
#   s.range <- range(pca$workGrid)
#   tt <- sort(unique(c(t)))
#   in.idx <- tt >= s.range[1] & tt <= s.range[2]
#   tt <- tt[in.idx]; m <- length(tt)
#   impute_res <- lapply(1:n, function(j) {
#     ConvertSupport(pca$workGrid, toGrid = tt, mu = pred[j,])
#   })
#   betat <- foreach(j = 1:m, .combine = "rbind") %dopar% {
#     tj <- tt[j]
#     xj <- sapply(impute_res, '[', j)
#     yj <- sapply(1:n, function(x) {
#       if (tj %in% t[x,]) {
#         Ly[[x]][which(t[x,]==tj)]
#         } else {
#         NA
#         }
#     })
#     nna.idx <- which(!is.na(yj))
#     yj <- yj[nna.idx]
#     Xj <- cbind(1, xj[nna.idx])
#     bj <- tryCatch(solve(t(Xj) %*% Xj) %*% t(Xj) %*% yj, error = function(e) {rep(NA, 2)})
#     c(bj)
#   }
#   return(list(beta = betat, time = tt, pca.sigma2 = pca$sigma2, pca.rho = pca$rho, pca.K = pca$selectK, 
#               pca.mubw = pca$bwMu, pca.covbw = pca$bwCov))
# }
# 
# ## function for last observation carried forward
# tdLOCFlassoX <- function(Ly, t, Ws, s) {
#   n <- nrow(t)
#   m <- ncol(t)
#   dx <- ncol(Ws)
#   nXhat <- matrix(0, ncol = dx, nrow = m*n)
#   for (v in 1:dx) {
#     Wsv <- matrix(Ws[,v], ncol = m, nrow = n, byrow = TRUE)
#     sv <- s[[v]]
#     nXhat[,v] <- c(sapply(1:n, function(i) {
#       sapply(t[i,], function(y) {last_obs(Wsv[i,], sv[i,], y, nearest = FALSE)})
#     }))
#   }
#   tt <- sort(unique(c(t)))
#   betat <- foreach(j = 1:length(tt), .combine = "rbind") %dopar% {
#     tj = tt[j]
#     yj <- unlist(lapply(1:n, function(x) {
#       if (tj %in% t[x,]) {
#         Ly[[x]][which(t[x,]==tj)]
#       }
#     }))
#     xj <- lapply(1:n, function(x) {
#       if (tj %in% t[x,]) {
#         nXhat[which(t[x,]==tj)+m*(x-1),,drop=FALSE]
#       }
#     })
#     xj <- do.call("rbind", xj)
#     nna.idx <- which(apply(xj, 1, function(x) {all(!is.na(x))}))
#     Xj <- cbind(1, xj[nna.idx,,drop=FALSE])
#     yj <- yj[nna.idx]
#     #bj <- solve(t(Xj) %*% Xj) %*% t(Xj) %*% yj
#     bj <- tryCatch(solve(t(Xj) %*% Xj) %*% t(Xj) %*% yj, error = function(e) {rep(NA, 2)})
#     c(bj) 
#   }
#   return(list(beta = betat, time = tt))
# }
# 
