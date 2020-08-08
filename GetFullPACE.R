## functions for getting FullPACE xi 

## univariate FPCA and substract estimated mean from observation
uFPCACov <- function(W_s, s, n, m, optns) {
  dx <- ncol(W_s)
  pca.list = foreach(v = 1:dx) %dopar% {
    # do FPCA, impute W_s on T domain 
    Ly <- split(W_s[, v], as.factor(rep(1:n, each = m)))
    Lt <- lapply(1:n, function(x) {s[[v]][x,]})
    pca <- FPCAIC(Ly, Lt, optns, spec = list("AIC", FALSE))
    pca
  }
  ## extract resid list
  resid.list <- lapply(pca.list, "[[", "resid")
  ## reformat W_i into a list of length n: each is ordered by variable
  W_is <- lapply(1:n, function(x) {
    unlist(lapply(1:dx, function(v) {resid.list[[v]][[x]]}))
  })
  ## cross covariance estimation
  cr.cov <- foreach(v1 = 1:(dx-1)) %dopar% {
    cov.v1 <- list()
    for (v2 in (v1+1):dx) {
      Ly1 <- resid.list[[v1]]
      Ly2 <- resid.list[[v2]]
      Lt1 <- lapply(1:n, function(x) {s[[v1]][x,]})
      Lt2 <- lapply(1:n, function(x) {s[[v2]][x,]})
      cov <- GetCrCovYX(Ly1 = Ly1, Lt1 = Lt1, Ly2 = Ly2, Lt2 = Lt2, Ymu1 = rep(0, length(unique(unlist(Lt1)))), 
                        Ymu2 = rep(0, length(unique(unlist(Lt2)))), bw1 = 1, bw2 = 1, kern = "epan",
                        grid1 = pca.list[[v1]]$workGrid, grid2 = pca.list[[v2]]$workGrid)
      cov.v1[[v2]] <- cov
    }
    cov.v1
  }
  return(list(pca.list = pca.list, resid.list = resid.list, cr.cov = cr.cov, W_is = W_is))
}

## estimate cross omegas, currently, fix K = 3 for all variables
GetOmega_v <- function(cr.cov, pca.list) {
  dx <- length(pca.list)
  Omegas <- list()
  for (v1 in 1:(dx-1)) {
    omega.v1 <- list()
    for (v2 in (v1+1):dx) {
      cov.obj <- cr.cov[[v1]][[v2]]
      cov <- cov.obj$smoothedCC # column is corresponding to v1 and row is corresponding to v2
      # grid <- cov.obj$smoothGrid # the first column is v1, the second is v2
      # eigenfunctions: should be in the same grid as the grid
      phi1 <- pca.list[[v1]]$phi; grid1 <- pca.list[[v1]]$workGrid
      phi2 <- pca.list[[v2]]$phi; grid2 <- pca.list[[v2]]$workGrid
      # in case there are numerical issues: impute eigenfunction onto covariance matrix
      # grid1[1] <- grid[1,1]
      # grid2[1] <- grid[1,2]
      # grid1[grid.num] <- grid[grid.num,1]
      # grid2[grid.num] <- grid[grid.num,2]
      # phi1 <- ConvertSupport(fromGrid = grid1, toGrid = grid[,1], phi = phi1)
      # phi2 <- ConvertSupport(fromGrid = grid2, toGrid = grid[,2], phi = phi2)
      omega12 <- matrix(0, ncol = ncol(phi2), nrow = ncol(phi1))
      for (k1 in 1:ncol(phi1)) {
        for (k2 in 1:ncol(phi2)) {
          temp <- unlist(sapply(1:length(grid2), function(x) {trapzRcpp(X = grid1, Y = phi1[,k1]*cov[,x])}))
          omega12[k1,k2] <- trapzRcpp(X = grid2, Y = phi2[,k2]*temp)
        }
      }
      omega.v1[[v2]] <- omega12
      ## update cross covariance with fitted covariance
      cr.cov[[v1]][[v2]][['fittedCC']] <- phi1%*%omega12%*%t(phi2)
    }
    Omegas[[v1]] <- omega.v1
  }
  # omega matrix for the same variable
  for (v in 1:(dx-1)) {
    Omegas[[v]][[v]] <- diag(c(pca.list[[v]]$lambda))
  }
  Omegas[[dx]] <- list()
  Omegas[[dx]][[dx]] <- diag(c(pca.list[[dx]]$lambda))
  ## reformat and construct the Omega matrix into an array
  Omega_vs <- foreach(v = 1:dx) %dopar% {
    Omegav <- NULL
    for (i in 1:dx) {
      if (i < v) {
        Omegav <- cbind(Omegav, t(Omegas[[i]][[v]]))
      } else {
        Omegav <- cbind(Omegav, Omegas[[v]][[i]])
      }
    }
    Omegav
  }
  return(list(Omega_vs = Omega_vs, Omegas = Omegas, cr.cov = cr.cov))
}

## function for constructing Phi
GetPhis <- function(pca.list, n = 100) {
  dx <- length(pca.list)
  Phis <- foreach(i = 1:n) %dopar% {
    phi.list <- lapply(pca.list, function(x) {t(x$phiObs[[i]])})
    as.matrix(bdiag(phi.list))
  }
  return(Phis)
}

## convert support for cross covariance
ConvertSuportcr <- function(fromGrid, toGrid, Cov = NULL) {
  for (i in 1:2) {
    fgrid <- fromGrid[,i]; tgrid <- toGrid[,i]
    # In case the range of toGrid is larger than fromGrid due to numeric error
    #buff <- .Machine$double.eps * max(abs(fromGrid)) * 10
    buff <- 10e-06
    if (abs(tgrid[1] - fgrid[1]) < buff)
      tgrid[1] <- fgrid[1]
    if (abs(tgrid[length(tgrid)] - fgrid[length(fgrid)]) < buff)
      tgrid[length(tgrid)] <- fgrid[length(fgrid)]
    if ( ( fgrid[1] - buff  >  tgrid[1]) || 
         ( fgrid[length(fgrid)] + buff < tgrid[length(tgrid)]) ) {
      stop("Insufficient size of 'fromGrid'.")}
    fromGrid[,i] <- fgrid; toGrid[,i] <- tgrid
  }
  gd <- expand.grid(X=toGrid[,1], Y=toGrid[,2])
  ret <- matrix(interp2lin(fromGrid[,1], fromGrid[,2], Cov, gd$X, gd$Y), nrow=length(toGrid[,1]))
  return(ret)
}
## function for constructing Sigma_Wi
GetSigma_Wi <- function(pca.list, cr.cov, s, n = 100, m = 5) {
  dx <- length(pca.list)
  Sigmas <- foreach(i = 1:n) %dopar% {
    # off diagonal matrices
    v1.list <- foreach(v1 = 1:(dx-1)) %do% {
      v2.list <- foreach(v2 = (v1+1):dx) %do% {
        cov <- cr.cov[[v1]][[v2]]$fittedCC
        grid <- cr.cov[[v1]][[v2]]$smoothGrid
        s1 <- s[[v1]][i,]
        s2 <- s[[v2]][i,]
        ConvertSuportcr(grid, cbind(s1, s2), cov)
        #gd <- expand.grid(x = s1, y = s2)
        #matrix(interp2lin(xin = grid[,1], yin = grid[,2], zin = cov, xou = gd$x, you = gd$y), nrow = m)
      }
      cbind(matrix(0, ncol = m*v1, nrow = m), do.call("cbind", v2.list))
    }
    off.mat <- rbind(do.call("rbind", v1.list), matrix(0, ncol = m*dx, nrow = m))
    mat <- (off.mat + t(off.mat))
    # digonal matrices
    for (v in 1:dx) {
      vidx <- (m*(v-1)+1):(m*v)
      mat[vidx, vidx] <- pca.list[[v]]$fittedCovObs[[i]] + diag(pca.list[[v]]$sigma2, m)
    }
    # for (v in 1:dx) {
    #   #cov <- pca.list[[v]]$fittedCov #+ diag(pca.list[[v]]$sigma2, grid.num)
    #   cov <- pca.list[[v]]$phi[,1:K] %*% diag(pca.list[[v]]$lambda[1:K]) %*% t(pca.list[[v]]$phi[,1:K])
    #   grid <- pca.list[[v]]$workGrid
    #   # in case there are numerical issues
    #   ss <- s[[v]][i,]
    #   if (ss[1] < min(grid)) {
    #     grid[1] <- ss[1]
    #   }
    #   if (ss[m] > max(grid)) {
    #     grid[grid.num] <- ss[m]
    #   }
    #   vidx <- (m*(v-1)+1):(m*v)
    #   gd <- expand.grid(x = ss, y = ss)
    #   mat[vidx, vidx] <- matrix(interp2lin(xin = grid, yin = grid, zin = cov, xou = gd$x, you = gd$y), nrow = length(ss)) + diag(pca.list[[v]]$sigma2, m)
    # }
    mat
  }
  # Sigmas_inv <- lapply(Sigmas, function(x) {
  #   eig <- eigen(x, symmetric = TRUE)
  #   pos.idx <- eig$values > 0
  #   eig$vectors[,pos.idx,drop=FALSE] %*% diag(1/eig$values[pos.idx]) %*% t(eig$vectors[,pos.idx,drop=FALSE])
  # })
  Sigmas_inv <- lapply(Sigmas, function(x) {
    solve(x)
  })
  return(Sigmas_inv)
}

## calculate pc score conditioning on all observations
GetFullPCScore <- function(pca.list, cr.cov, W_is, s, n = 100, m = 5) {
  Omega_obj <- GetOmega_v(cr.cov, pca.list)
  # update cross covariance
  cr.cov1 <- Omega_obj$cr.cov
  Omega_vs <- Omega_obj$Omega_vs
  Phis <- GetPhis(pca.list, n)
  Sigmas_inv <- GetSigma_Wi(pca.list, cr.cov1, s, n, m)
  ## estimate xi_iv: PC score for each subject and each variable
  dx <- length(pca.list)
  xi.list <- foreach(v = 1:dx) %dopar% {
    xi <- lapply(1:n, function(i) {
      t(Omega_vs[[v]] %*% Phis[[i]] %*% Sigmas_inv[[i]] %*% W_is[[i]])
    })
    do.call("rbind", xi)
  }
  return(xi.list)
} 


