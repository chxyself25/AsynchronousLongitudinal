## lasso regression 
lasso_2step <- function(MX, Y, type, pars.n, intercept = TRUE, nfolds) {
  # lasso regression
  if (type %in% c('mse', 'mae')) {
    fit.fm = cv.glmnet(x=MX, y=Y, family = "gaussian", type.measure = type, lambda = seq(0.01, 1, length.out = 100),
                       nfolds = nfolds, alpha = 1, intercept = intercept, standardize = FALSE)
    beta.fm <- unlist(coef(fit.fm, s = 'lambda.1se'))
  } else if (type == 'BIC'){ 
    fit.fm <- cv.glmnet(x=MX, y=Y, family = "gaussian", type.measure = 'mse', 
                        nfolds = nfolds, alpha = 1, standardize = FALSE, intercept = FALSE)
    bicss <- sapply(1:length(fit.fm$lambda), function(x) {
      k <- sum(unlist(coef(fit.fm, s = fit.fm$lambda[x])) != 0)
      sigma2 <- fit.fm$cvm[x]*(n*m)/(n*m-k)
      log(n*m)*k + n*m*log(2*pi*sigma2) + n*m-k
    })
    beta.fm <- unlist(coef(fit.fm, s = fit.fm$lambda[which.min(bicss)]))[-1]
  }
  c.idx <- which(beta.fm[-1] != 0) # correct variables identified except for intercept
  ic.idx <- which(beta.fm != 0)
  if (length(c.idx) == 0) { # all variables are zero
    betas <- rep(NA, pars.n+1)
  } else {
    betas <- rep(0, pars.n+1)
    ones <- ifelse(beta.fm[1] != 0, 1, NULL)
    nMX <- cbind(ones, MX[,c.idx])
    betas[ic.idx] <- solve(t(nMX) %*% nMX) %*% t(nMX) %*% Y
  }
  return(list(betas = betas, c.idx = c.idx))
}

## simplified 2 step lasso
lasso_2step0 <- function(MX, Y, pars.n, nfolds) {
  # lasso regression
  fit.fm = cv.glmnet(x=MX, y=Y, family = "gaussian", type.measure = "mse", lambda = seq(0.01, 1, length.out = 100),
                     nfolds = nfolds, alpha = 1, intercept = FALSE, standardize = FALSE)
  beta.fm <- unlist(coef(fit.fm, s = 'lambda.1se'))
  c.idx <- which(beta.fm[-1] != 0) # correct variables identified except for intercept
  #ic.idx <- which(beta.fm[-1] != 0)
  if (length(c.idx) == 0) { # all variables are zero, not included
    betas <- rep(0, pars.n)
  } else {
    betas <- rep(0, ncol(MX))
    nMX <- MX[,c.idx,drop=FALSE]
    betas[c.idx] <- solve(t(nMX) %*% nMX) %*% t(nMX) %*% Y
  }
  return(list(betas = betas, c.idx = c.idx))
}

# adaptive lasso for functional match: type is mse only currently
adapt_lasso <- function(MX, Y, intercept = FALSE, nfolds = 10) {
  if (!intercept) {
    b0 <- solve(t(MX)%*%MX)%*%t(MX)%*%Y
    fit.fm = cv.glmnet(x=MX, y=Y, family = "gaussian", type.measure = "mse", lambda = seq(0.001, 1, length.out = 100),
                       nfolds = nfolds, alpha = 1, intercept = FALSE, standardize = FALSE, penalty.factor = 1/abs(b0))
    beta.fm <- unlist(coef(fit.fm, s = 'lambda.1se'))
    c.idx <- which(beta.fm[-1] != 0)
  } else {
    b0 <- solve(t(cbind(1,MX))%*%cbind(1,MX))%*%t(cbind(1,MX))%*%Y
    fit.fm = cv.glmnet(x=MX, y=Y, family = "gaussian", type.measure = "mse", lambda = seq(0.001, 1, length.out = 100),
                       nfolds = nfolds, alpha = 1, intercept = TRUE, standardize = FALSE, penalty.factor = 1/abs(b0))
    beta.fm <- unlist(coef(fit.fm, s = 'lambda.1se'))
    c.idx <- which(beta.fm[-1] != 0)
  }
  return(list(betas = beta.fm, c.idx = c.idx))
}

