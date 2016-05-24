
KF <- function(y, ss, convergence = c(0.001, length(y)), t0 = 1)
{
  stopifnot(convergence[2] >= 1)
  
  n <- length(y)
  p <- 1 # ncol(y) univariate data

  a0 <- ss$a0
  P0 <- ss$P0
  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  tZ <- t(Z)
  tmT <- t(mT)

  notconv <- TRUE
  counter <- 0
  tol <- convergence[1]
  maxiter <- convergence[2]
  checkconv <- maxiter < length(y)
  convit <- if (checkconv) 0 else NULL

  # storage vectors

  a.pred <- matrix(nrow = n, ncol = length(a0))
  P.pred <- array(NA, dim = c(dim(P0), n))
  a.upd <- matrix(nrow = n + 1, ncol = length(a0))
  P.upd <- array(NA, dim = c(dim(P0), n + 1))
  K <- matrix(nrow = n, ncol = length(a0))
  L <- array(NA, dim = c(length(a0), length(a0), n))
  v <- rep(NA, n)
  f <- rep(NA, n)
  lik.contrib <- rep(0, n)

  a.upd[1,] <- a0
  P.upd[,,1] <- P0

  # filtering recursions

  for (i in seq_len(n))
  {
    a.pred[i,] <- mT %*% a.upd[i,]
    
    if (notconv) {
      P.pred[,,i] <- mT %*% P.upd[,,i] %*% tmT + Q
    } else 
      P.pred[,,i] <- P.pred[,,i-1]

    if (!is.na(y[i]))
    {
      # prediction

      v[i] <- y[i] - Z %*% a.pred[i,] # prediction error
      
      if (notconv) {
        Mt <- tcrossprod(P.pred[,,i], Z)

        f[i] <- Z %*% Mt + H # variance of prediction error ('gain' in arima.c)
      } else 
        f[i] <- f[i-1]

      # contribution to the log-likelihood

      lik.contrib[i] <- log(f[i]) + v[i]^2 / f[i]

      # updating

      K[i,] <- Mt / f[i]

      a.upd[i+1,] <- a.pred[i,] + K[i,] * v[i]

      if (notconv) {
        P.upd[,,i+1] <- P.pred[,,i] - tcrossprod(Mt) / f[i]

        # required in this form by disturbance smoother
        K[i,] <- mT %*% K[i,]
        # required for Kalman smoother
        L[,,i] <- mT - K[i,] %*% Z
      } else {
        P.upd[,,i+1] <- P.upd[,,i]
        K[i,] <- K[i-1,]
        L[,,i] <- L[,,i-1]
      }

      # check convergence of the filter:
      if (checkconv && notconv)
      {
        if (i == 1) {
          fprev <- f[i] + tol + 1
        }
      
        if (abs(f[i] - fprev) < tol)
        {
          # remain steady over 'maxiter' consecutive iterations
          if (convit == i - 1)
          {
            counter <- counter + 1
          } else
            counter <- 1
          convit <- i
        }

        fprev <- f[i]

##FIXME this should be insice the if statement where counter is incremented
        if (counter == maxiter) {
          notconv <- FALSE # the filter has converged to a steady state
          convit <- i
        }
      }
    } else # if (is.na(y[i]))
    {
      a.upd[i+1,] <- a.pred[i,]
      P.upd[,,i+1] <- P.pred[,,i]
      #v[i] <- NA
      #f[i] <- NA #P.pred + H

      # NOTE reset if NA is found after the filter converged
      # for safety and also because it is a way to deal with the NA in f[i]
      if (!notconv) {
        notconv <- TRUE
        counter <- 1
      }
    }
  }

  if (notconv)
    convit <- NULL

  v <- ts(v)
  f <- ts(f)
  a.upd <- ts(a.upd[-1,])
  K <- ts(K)
  tsp(v) <- tsp(f) <- tsp(a.upd) <- tsp(K) <- tsp(y)

  mll <- 0.5 * (n - t0 + 1) * log(2 * pi) + 
    0.5 * sum(lik.contrib[seq.int(t0, n)], na.rm = TRUE)

  list(v = v, f = f, K = K, L = L, a.upd = a.upd, P.upd = P.upd[,,-1],
    a.pred = a.pred, P.pred = P.pred, mll = mll, convit = convit)
}

KF.C <- function(y, ss, convergence = c(0.001, length(y)), t0 = 1) #sUP = 1
{
  stopifnot(convergence[2] >= 1)
  checkconv <- as.integer(convergence[2] < length(y))

  y[is.na(y)] <- -9999.99

  #NOTE non-diagonal matrices are passed transposed
  mll <- .C("KF_C", 
    dim = as.integer(c(length(y), ncol(ss$Z), ncol(ss$V), t0, checkconv)), #sUP
    y = y, sZ = as.numeric(ss$Z), sT = as.numeric(t(ss$T)), 
    H = as.numeric(ss$H), #sV = as.numeric(ss$V), 
    sQ = as.numeric(ss$Q), sa0 = as.numeric(ss$a0), sP0 = as.numeric(ss$P0),
    convtol = as.numeric(convergence[1]), convmaxiter = as.integer(convergence[2]),
    mll = double(1), PACKAGE = "KFKSDS")$mll

  mll
}

KF.deriv <- function(y, ss, xreg = NULL, convergence = c(0.001, length(y)), t0 = 1)
{
  n <- length(y)
  r <- ncol(ss$V)
  rp1 <- r + 1

  a0 <- ss$a0
  P0 <- ss$P0
  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  tZ <- t(Z)
  tmT <- t(mT)

  notconv <- TRUE
  counter <- 0
  tol <- convergence[1]
  maxiter <- convergence[2]
  checkconv <- maxiter < n
  convit <- if (checkconv) 0 else NULL

  if (!is.null(xreg))
  {
    # no check about the correct definition of "xreg"
    # which should be a list containg a matrix (xreg) with 
    # length(y) number of rows and a vector with the 
    # corresponding coefficients (coefs)
    if (is.list(xreg)) {
      y <- y - xreg$xreg %*% cbind(xreg$coefs)
      xreg <- xreg$xreg
    } # otherwise "y" is assumed to be passed already as y-xreg*coefs
    ncxreg <- ncol(xreg)    
  } else
    ncxreg <- 0

  # storage vectors

  a.pred <- matrix(nrow = n, ncol = length(a0))
  P.pred <- array(NA, dim = c(dim(P0), n))
  a.upd <- matrix(nrow = n + 1, ncol = length(a0))
  P.upd <- array(NA, dim = c(dim(P0), n + 1))

  varnms <- paste("var", seq_len(rp1), sep = "")
  da.pred <- array(dim = c(n, ncol(Z), rp1 + ncxreg))
  dimnames(da.pred)[[3]] <- c(varnms, colnames(xreg))
  dP.pred <- array(dim = c(dim(P0), n, rp1))
  dimnames(dP.pred)[[4]] <- varnms
  da.upd <- array(0, dim = c(n + 1, ncol(Z), rp1 + ncxreg))
  dimnames(da.upd)[[3]] <- dimnames(da.pred)[[3]]
  dP.upd <- array(0, dim = c(dim(P0), n + 1, rp1))
  dimnames(dP.upd)[[4]] <- varnms
  dv <- matrix(nrow = n, ncol = rp1 + ncxreg)
  colnames(dv) <- dimnames(da.pred)[[3]]
  df <- matrix(nrow = n, ncol = rp1, dimnames = list(NULL, varnms))
  dK <- array(dim = c(n, ncol(Z), rp1), dimnames = list(NULL, NULL, varnms))

  K <- matrix(nrow = n, ncol = length(a0))
  L <- array(NA, dim = c(length(a0), length(a0), n))
  v <- rep(NA, n)
  f <- rep(NA, n)
  lik.contrib <- rep(0, n)

  a.upd[1,] <- a0
  P.upd[,,1] <- P0

  #for (i in seq_len(rp1))
  #{
  #  da.upd[1,,i] <- 0
  #  dP.upd[,,1,i] <- 0
  #}

  # filtering recursions

  for (i in seq_len(n))
  {
    ip1 <- i + 1

    a.pred[i,] <- mT %*% a.upd[i,]
    P.pred[,,i] <- mT %*% P.upd[,,i] %*% tmT + Q

    if (!is.na(y[i]))
    {
      # prediction

      v[i] <- y[i] - Z %*% a.pred[i,] # prediction error
    
      if (notconv) {
        Mt <- tcrossprod(P.pred[,,i], Z)
        f[i] <- Z %*% Mt + H
      } else 
        f[i] <- f[i-1]
    
      # contribution to the log-likelihood

      lik.contrib[i] <- log(f[i]) + v[i]^2 / f[i]

      # updating

      K[i,] <- Mt / f[i]

      a.upd[ip1,] <- a.pred[i,] + K[i,] * v[i]
    
      if (notconv) {
        P.upd[,,ip1] <- P.pred[,,i] - tcrossprod(Mt) / f[i]
        K[i,] <- mT %*% K[i,]
        # required for Kalman smoother
        L[,,i] <- mT - K[i,] %*% Z
      } else {
        P.upd[,,ip1] <- P.upd[,,i]
        K[i,] <- K[i-1,]
        L[,,i] <- L[,,i-1]
      }

      # check convergence of the filter

      if (checkconv && notconv)
      {
        # NA not allowed in the first observation
        if (i == 1) {
          fprev <- f[i] + tol + 1
        }

        ref <- abs(f[i] - fprev)
        if (is.na(ref))
          ref <- tol + 1
        if (ref < tol)
        {
          # remain steady over 'maxiter' consecutive iterations
          if (convit == i - 1)
          {
            counter <- counter + 1
          } else
            counter <- 1
          convit <- i
        }

        fprev <- f[i]

        if (counter == maxiter) {
          notconv <- FALSE # the filter has converged to a steady state
          convit <- i
        }
      }
    } else { # if (is.na(y[i]))    
      a.upd[ip1,] <- a.pred[i,]
      P.upd[,,ip1] <- P.pred[,,i]
      #v[i] <- NA
      #f[i] <- NA #P.pred + H
    }

    # derivative terms

    fsq <- f[i]^2

    for (j in seq_len(rp1))
    {
      da.pred[i,,j] <- mT %*% da.upd[i,,j]
      dv[i,j] <- -Z %*% da.pred[i,,j]

      if (notconv)
      {
        if (j == 1) {
          dP.pred[,,i,j] <- mT %*% dP.upd[,,i,j] %*% tmT
        } else {
          if (r == 1) {
            tmp <- matrix(1)
          } else {
            tmp <- diag(rep(0, r))
            diag(tmp)[j-1] <- 1
          }
          tmp <- ss$R %*% tmp %*% t(ss$R)
          dP.pred[,,i,j] <- mT %*% dP.upd[,,i,j] %*% tmT + tmp
        }

        if (j == 1) {
          df[i,j] <- Z %*% dP.pred[,,i,j] %*% tZ + 1
        } else
          df[i,j] <- Z %*% dP.pred[,,i,j] %*% tZ

        dK[i,,j] <- mT %*% dP.pred[,,i,j] %*% tZ / f[i] - mT %*% P.pred[,,i] %*% tZ * df[i,j] / fsq

        dP.upd[,,ip1,j] <- dP.pred[,,i,j] - dP.pred[,,i,j] %*% tZ %*% t(Mt) / f[i] +
          tcrossprod(Mt) * df[i,j] / fsq - Mt %*% Z %*% dP.pred[,,i,j] / f[i]
      } else {
        im1 <- i - 1
        dP.pred[,,i,j] <- dP.pred[,,im1,j]
        df[i,j] <- df[im1,j]
        dK[i,,j] <- dK[im1,,j]
        dP.upd[,,ip1,j] <- dP.upd[,,i,j]
      }

      da.upd[ip1,,j] <- da.pred[i,,j] + dP.pred[,,i,j] %*% tZ * v[i] / f[i] -
        Mt * df[i,j] * v[i] / fsq + Mt * dv[i,j] / f[i]

      if (is.na(y[i])) {
        dK[i,,j] <- 0
        da.upd[ip1,,j] <- da.pred[i,,j]
        dP.upd[,,ip1,j] <- dP.pred[,,i,j]
      }
    } # end for (j in seq_len(rp1))

    if (!is.null(xreg))
    {
      k <- 1
      for (j in seq(rp1 + 1, ncol(dv)))
      {
        da.pred[i,,j] <- mT %*% da.upd[i,,j]
        dv[i,j] <- -Z %*% da.pred[i,,j] - xreg[i,k]
        da.upd[ip1,,j] <- da.pred[i,,j] + Mt * dv[i,j] / f[i]
        k <- k + 1

        if (is.na(y[i])) {
          da.upd[ip1,,j] <- da.pred[i,,j]
        }
      }
    }

    if (is.na(y[i]))
    {
      if (!notconv) {
        notconv <- TRUE
        counter <- 1
      }

      dv[i,] <- 0
      df[i,] <- 0
    }
  }

  if (notconv)
    convit <- NULL

  v <- ts(v)
  f <- ts(f)
  a.upd <- ts(a.upd[-1,])
  K <- ts(K)
  tsp(v) <- tsp(f) <- tsp(a.upd) <- tsp(K) <- tsp(y)

  dvof <- (dv[,varnms]*c(f) - c(v)*df) / c(f)^2

  m <- ncol(ss$Z)
  dL <- array(0, dim = c(m, m, rp1, n))
  dimnames(dL)[[3]] <- varnms #dimnames(da.pred)[[3]]
  for (k in seq_len(dim(dL)[3]))
  {
    for (i in seq_len(n))
    { 
      dL[,,k,i] <- - dK[i,,k] %*% Z
    }
  }

  mll <- 0.5 * (n - t0 + 1) * log(2 * pi) + 
    0.5 * sum(lik.contrib[seq.int(t0, n)], na.rm = TRUE)

  list(v = v, f = f, K = K, L = L, a.upd = a.upd, P.upd = P.upd[,,-1],
    a.pred = a.pred, P.pred = P.pred, 
    da.pred = da.pred, dP.pred = dP.pred,
    da.upd = da.upd[-1,,], dP.upd = dP.upd[,,-1,],
    dv = dv, #dv = colSums(dv), 
    df = df, #df = colSums(df), 
    dvof = dvof, 
    dK = dK, #dK = colSums(dK), 
    dL = dL, mll = mll, convit = convit)
}

KF.deriv.C <- function(y, ss, xreg = NULL, convergence = c(0.001, length(y)), 
  t0 = 1, return.all = FALSE)
{
  stopifnot(convergence[2] >= 1)
  checkconv <- as.integer(convergence[2] < length(y))

  n <- length(y)
  r <- ncol(ss$V)
  m <- ncol(ss$Z)
  nm <- n * m
  nmm <- nm*m
  nr <- n * r
  rp1 <- r + 1
  nrp1 <- nr + n #n * rp1
  nrp1m <- nrp1 * m

  if (!is.null(xreg))
  {
    # no check about the correct definition of "xreg"
    # which should be a list containg a matrix (xreg) with 
    # length(y) number of rows and a vector with the 
    # corresponding coefficients (coefs)
    if (is.list(xreg)) {
      y <- y - xreg$xreg %*% cbind(xreg$coefs)
      xreg <- xreg$xreg
    } # otherwise "y" is assumed to be passed already as y-xreg*coefs
    ncxreg <- ncol(xreg)    
  } else
    ncxreg <- 0

  y[is.na(y)] <- -9999.99

  res <- .C("KF_deriv_C", 
    dim = as.integer(c(n, m, r, t0, checkconv, ncxreg)), #sUP
    y = as.numeric(y), xreg = as.numeric(xreg), 
    sZ = as.numeric(ss$Z), sT = as.numeric(t(ss$T)), 
    sH = as.numeric(ss$H), sR = as.numeric(t(ss$R)), sV = as.numeric(ss$V), 
    sQ = as.numeric(ss$Q), sa0 = as.numeric(ss$a0), sP0 = as.numeric(ss$P0),
    convtol = as.numeric(convergence[1]), convmaxiter = as.integer(convergence[2]),
    conv = integer(2), mll = double(1), invf = double(n), vof = double(n),
    dvof = double(nrp1), dfinvfsq = double(nrp1),
    dv = double(nrp1 + n*ncxreg), df = double(nrp1), 
    a_pred0 = double(nm), P_pred0 = double(nmm),
    K0 = double(nm), L0 = double(nmm), da_pred0 = double(nrp1m + nm*ncxreg),
    dP_pred0 = double(nrp1m*m), dK0 = double(nrp1m),
    PACKAGE = "KFKSDS")

  df <- matrix(res$df, nrow = n, ncol = rp1)
  colnames(df) <- paste("var", seq_len(ncol(df)), sep = "")
  dv <- matrix(res$dv, nrow = n, ncol = rp1 + ncxreg)
  colnames(dv) <- c(colnames(df), colnames(xreg))

  if (res$conv[1] == 0) {
    convit <- res$conv[2] + 1
  } else 
    convit <- NULL

  # required for Kalman smoother

  if (return.all)
  {
    L <- array(res$L, dim = c(m, m, n))
    dK0 <- res$dK0  
    da.pred <- array(dim = c(n, m, ncol(dv)))
    dK <- array(dim = c(n, m, rp1))
    rp1m <- rp1 * m
    ref1 <- seq_len(rp1m)
    dP.pred <- array(dim = c(m, m, n, rp1))
    dPpred0 <- res$dP_pred0
    ref2 <- seq_len(m*m)
    for (i in seq_len(n))
    {
      L[,,i] <- t(L[,,i])
      dK[i,,] <- dK0[ref1]
      dK0 <- dK0[-ref1]

      for (j in seq_len(rp1))
      {
        # it is a symmetric matrix, byrow = TRUE or FALSE are possible
        dP.pred[,,i,j] <- matrix(dPpred0[ref2], nrow = m, ncol = m)
        dPpred0 <- dPpred0[-ref2]
      }
    }

    dapred0 <- res$da_pred0
    ref <- seq_len(nm)
    for (i in seq.int(rp1 + ncxreg))
    {
      da.pred[,,i] <- matrix(dapred0[ref], nrow = n, ncol = m, byrow = TRUE)
      dapred0 <- dapred0[-ref]
    }

    dL <- array(0, dim = c(m, m, rp1, n))
    dimnames(dL)[[3]] <- dimnames(da.pred)[[3]]
    for (k in seq(dim(dL)[3]))
    {
      for (i in seq(along = y))
      { 
        dL[,,k,i] <- - dK[i,,k] %*% ss$Z
      }
    }
  } else
    L <- dL <- da.pred <- dP.pred <- dK <- NULL

  f <- ts(1 / res$invf)
  v <- ts(res$vof * f)
  #dvof <- (dv*c(f) - c(v)*df) / c(f)^2
  dvof <- matrix(res$dvof, nrow = n, ncol = rp1)
  a.pred <- ts(matrix(res$a_pred, nrow = n, ncol = m, byrow = TRUE))
  K <- ts(matrix(res$K, nrow = n, ncol = m, byrow = TRUE))
  tsp(v) <- tsp(f) <- tsp(a.pred) <- tsp(K) <- tsp(y)

  list(a.pred = a.pred, P.pred = array(res$P_pred, dim = c(m, m, n)),
    v = v, f = f, K = K, L = L, da.pred = da.pred, dP.pred = dP.pred, 
    dv = dv, df = df, dvof = dvof, #invf = res$invf,
    dK = dK, dL = dL, mll = res$mll, convit = convit)
}
