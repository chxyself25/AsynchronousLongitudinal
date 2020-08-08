# function that can return the nearest observation index of s
# no last observation carred forward at all
nearest_obs <- function(s, tt) {
  return(which.min(abs(s-tt)))
}

## functions that may not be used often this moment
# function for calculating sigma0 estimates
# pca.list is pca results for each variable, 
# t is y observed points, min.t and max.t constructs range of t
# s is a list of x observations for all variables, has the same range with t
sigma0hat <- function(W_s, pca.list, t, min.t, max.t, s) {
  dx <- length(pca.list) # variables dimension
  n <- nrow(t) # number of subject
  m <- ncol(t) # number of observations for each subject, assume same
  K <- max(sapply(pca.list, function(x) {length(x$lambda)})) # maximum number of principle components as constant K
  bk1 <- m*K # level on dimension: dimension in each block of Li
  # calculate all the cross-covariance: omega_jl
  acomb3 <- function(...) abind::abind(..., along=3)
  Omega_jl <- foreach(j = 1:(dx-2), .multicombine=TRUE) %dopar% {
    foreach(l = (j+1):dx, .combine = 'acomb3', .multicombine = TRUE) %dopar% {
      pcaj <- pca.list[[j]]
      pcal <- pca.list[[l]]
      omegajl <- omega_jl(W_s=W_s, s=s, j=j, l=l, phij=pcaj$phi, phil=pcal$phi, tj=pcaj$workGrid, tl=pcal$workGrid,
                          m=m, n=n, K=K)
      return(omegajl)
    }
  }
  pcaj <- pca.list[[(dx-1)]]
  pcal <- pca.list[[dx]]
  Omega_jl[[(dx-1)]] <- abind(omega_jl(W_s=W_s, s=s, j=dx-1, l=dx, phij=pcaj$phi, phil=pcal$phi, tj=pcaj$workGrid, tl=pcal$workGrid,
                                       m=m, n=n, K=K), along=3)
  # go through each subject
  # sigmas <- 0
  sigmas = foreach(i = 1:n, .combine = 'acomb3') %dopar% {
    tt <- t[i,]
    Phii <- NULL
    Li <- matrix(0, ncol = m*K*dx, nrow = m*K*dx)
    temp.res <- list()
    for (v in 1:dx) {
      # convert phi to t domain and s domain
      ss <- s[[v]][i,]
      pcav <- pca.list[[v]]
      num.grid <- length(pcav$workGrid)
      phitv <- rbind(pcav$phi[1,], pcav$phi, pcav$phi[num.grid,])
      Phitv <- ConvertSupport(c(min.t, pcav$workGrid, max.t), toGrid = tt, phi = phitv)
      Phisv <- ConvertSupport(c(min.t, pcav$workGrid, max.t), toGrid = ss, phi= phitv)
      Lambdav <- diag(pcav$lambda)
      if (ncol(phitv) < K) {
        ghost0 <- matrix(0, ncol = K-ncol(Phitv), nrow = m)
        Phitv <- cbind(Phitv, ghost0)
        Phisv <- cbind(Phisv, ghost0)
        Lambdav <- diag(append(pcav$lambda, rep(0, K-ncol(phitv))))
      }
      # Sigma_{W_iv} on s domain
      sigma.iv <- Phisv %*% Lambdav %*% t(Phisv) + diag(rep(pcav$sigma2, m))
      Hiv <- Lambdav %*% t(Phisv) %*% solve(sigma.iv) %*% Phisv %*% Lambdav
      Li[(bk1*(v-1)+1):(bk1*v), (bk1*(v-1)+1):(bk1*v)] <- kronecker(Hiv, diag(rep(1,m)))
      Phii <- cbind(Phii, c(Phitv))
      temp.res[[v]] <- Lambdav %*% t(Phisv) %*% solve(sigma.iv) %*% Phisv
    }
    # calculate off diagonal elements of Li matrix
    for (d1 in 1:dx) {
      for (d2 in 1:dx) {
        if (d1 == d2) {
          next
        } else {
          j <- min(d1, d2)
          l <- max(d1, d2)
          Li[(bk1*(d1-1)+1):(bk1*d1), (bk1*(d2-1)+1):(bk1*d2)] <- kronecker(temp.res[[d1]]%*%Omega_jl[[j]][,,(l-j)], diag(rep(1, m)))
        }
      }
    }
    res <- c()
    for (j in 1:dx) {
      for (l in 1:dx) {
        res <- c(res, t(Phii[,j]) %*% Li[(bk1*(j-1)+1):(bk1*j), (bk1*(l-1)+1):(bk1*l)] %*% Phii[,l])
      }
    }
    matrix(res, byrow = TRUE, ncol = dx, nrow = dx)
  }
  return(apply(sigmas, 1:2, sum)/(n*m))
}
# function for calculating cross covariance, j =\= l
# W_s is observed x of all variables, each variable is one column
# s is a list of x observed points
# phij and phil are from pca results of each variable
# K is maximum number of principle components
omega_jl <- function(W_s, s, j, l, phij, phil, tj, tl, m, n, K) {
  Lyj <- split(W_s[,j], as.factor(rep(1:n, each = m)))
  Lyl <- split(W_s[,l], as.factor(rep(1:n, each = m)))
  Ltj <- lapply(1:n, function(x) {s[[j]][x,]})
  Ltl <- lapply(1:n, function(x) {s[[l]][x,]})
  num.obs1 <- length(unique(unlist(Ltj)))
  num.obs2 <- length(unique(unlist(Ltl)))
  smooth.cov <- GetCrCovYX(Ly1 = Lyj, Lt1 = Ltj, Ly2 = Lyl, Lt2 = Ltl, Ymu1 = rep(0, num.obs1), Ymu2 = rep(0, num.obs2), 
                           bw1=1, bw2=1)
  covjl <- smooth.cov$smoothedCC
  gridjl <- smooth.cov$smoothGrid # first column is for j and second is for l
  # numerical integration to get omega_jl
  omegajl <- matrix(0, ncol = K, nrow = K)
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      if (k1 > ncol(phij) | k2 > ncol(phil)) {
        omegajl[k1, k2] <- 0
      } else {
        phij <- ConvertSupport(fromGrid = c(0,tj,10), toGrid = gridjl[,1], phi = rbind(phij[1,], phij, phij[51,]))
        phil <- ConvertSupport(fromGrid = c(0,tl,10), toGrid = gridjl[,2], phi = rbind(phil[1,], phil, phil[51,]))
        temp <- unlist(sapply(1:nrow(gridjl), function(x) {trapzRcpp(X = gridjl[,1], Y = phij[,k1]*covjl[,x])}))
        omegajl[k1, k2] <- trapzRcpp(X = gridjl[,2], Y = phil[,k2]*temp)
      }
    }
  }
  return(omegajl)
}

## may be abandoned
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

## calculate weights based on minimum different between X time and Y time
obs_weights <- function(s, t) {
  dx <- length(s)
  n <- nrow(t)
  weights <- NULL
  for (v in 1:dx) {
    sv <- s[[v]]
    weightv <- foreach(i = 1:n, .combine = 'c') %dopar% {
      (sapply(t[i,], function(x) {min(abs(sv[i,]-x))}))
    }
    weights <- cbind(weights, weightv)
  }
  return(rowSums(weights))
}

## CFM: corrected functional match
CFMlasso <- function(W_s, t, s, pca.list, Y, Xhat) {
  n <- nrow(t)
  m <- ncol(t)
  dx <- ncol(W_s)
  Sigma0hat <- sigma0hat(W_s, pca.list, t, min.t=0, max.t=10, s)
  # Dantzig selector
  ds.cfm <- slim(Xhat, Y, method = "dantzig", lambda = seq(0.001, 0.5, length.out = 100), Sigma0 = Sigma0hat)
  # select estimation by BIC
  bicss <- sapply(1:length(ds.cfm$lambda), function(x) {
    k <- sum(ds.cfm$beta0[[x]] != 0)
    log(n*m)*k + sum((Y- Xhat%*%(ds.cfm$beta0[[x]]))^2)
  })
  beta.cfm <- ds.cfm$beta0[[which.min(bicss)]]
  c.idx <- which(beta.cfm != 0)
  betas <- rep(0, dx)
  betas[c.idx] <- solve(t(Xhat[,c.idx]) %*% Xhat[,c.idx]) %*% t(Xhat[,c.idx]) %*% Y
}

## MFPCA method: Multivariate Functional Match: MFM
MFMlasso <- function(W_s, s, t, Z, Y, type = "mse", nfolds = 10, pca.list = NULL) {
  m <- ncol(t)
  n <- nrow(t)
  dx <- ncol(W_s)
  xiZ <- NULL
  if (is.null(pca.list)) {
    pca.list <- list()
    for (v in 1:dx) {
      Ly <- split(W_s[, v], as.factor(rep(1:n, each = 5)))
      Lt <- lapply(1:n, function(x) {s[[v]][x,]})
      #userMu <- list(t = seq(0, 10, length.out = 100), mu = rep(0, 100))
      pca <- FPCA(Ly, Lt, optns = list(dataType = "Sparse", userBwMu = 1, FVEthreshold = 0.99))
      pca.list[[v]] <- pca
      xiZ <- cbind(xiZ, pca$xiEst)
    }
  } else {
    for(v in 1:dx) {xiZ <- cbind(xiZ, pca.list[[v]]$xiEst)}
  }
  ## M1,..., Mp
  Mp <- cumsum(sapply(pca.list, function(x) {x$selectK}))
  Zhat <- t(xiZ)%*%xiZ/(n-1)
  zeig <- eigen(Zhat)
  pos.idx <- zeig[['values']] >= 0
  vhat <- zeig$values[pos.idx]
  chat <- zeig$vectors[, pos.idx, drop = FALSE]
  ## estimate eigenfunctions on the workGrid from FPCA
  mphi <- foreach(v = 1:dx) %do% {
    phiv <- pca.list[[v]]$phi
    if (v == 1) {chatv <- chat[1:(Mp[v]), ]} else{chatv <- chat[(Mp[v-1]+1):(Mp[v]), ]}
    phiv %*% chatv
  }
  ## principe component scores
  rhohat <- xiZ %*% chat
  ## imputation to t domain
  Xhat <- foreach(i = 1:n, .combine = 'rbind') %dopar% {
    tt <- t[i,]
    # # mean function
    # mu_t <- sapply(pca.list, function(x) {ConvertSupport(fromGrid = c(0,x$workGrid,10), toGrid = tt, 
    #                                              mu = c(x$mu[1], x$mu, x$mu[51]))})
    # # eigenfunctions
    # phi_t <- lapply(1:dx, function(v) {ConvertSupport(fromGrid = c(0, pca.list[[v]]$workGrid, 10), toGrid = tt,
    #                                                   phi = rbind(mphi[[v]][1,], mphi[[v]], mphi[[v]][51,]))})
    # mu_t + sapply(phi_t, function(x) {x %*% rhohat[i,]})
    impute_res <- lapply(1:dx, function(v) {
      pca.list[[v]]$mu + mphi[[v]] %*% rhohat[i,]
    })
    sapply(1:dx, function(v) {ConvertSupport(fromGrid = c(0, pca.list[[v]]$workGrid, 10), toGrid = tt,
                                             mu = c(impute_res[[v]][1], impute_res[[v]], impute_res[[v]][51]))})
  }
  MX <- cbind(Xhat, Z)
  if(is.null(Z)) {pars.n <- dx}else {pars.n <- dx + ncol(Z)}
  lasso.res <- lasso_2step(MX, Y, type = type, pars.n, intercept = TRUE, nfolds)
  return(list(Xhat=Xhat, pca.list = pca.list, betas = lasso.res$betas, selection = lasso.res$c.idx))
}
