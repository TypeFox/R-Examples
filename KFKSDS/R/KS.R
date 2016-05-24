
KS <- function(y, ss, kf)
{
  n <- length(y)
  if (is.null(kf$convit)) {
    convit <- n + 1
    nmconvit <- n # no effect inthe if statement below
  } else {
    convit <- kf$convit
    nmconvit <- n - convit
  }

  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  tZ <- t(Z)
  tmT <- t(mT)

  v <- kf$v
  a.pred <- kf$a.pred
  P.pred <- kf$P.pred
  f <- kf$f
  L <- kf$L

  r <- matrix(nrow = n + 1, ncol = ncol(a.pred))
  N <- array(dim = c(ncol(a.pred), ncol(a.pred), n + 1))
  ahat <- matrix(nrow = n, ncol = ncol(a.pred))
  varahat <- array(NA, dim = c(rep(ncol(a.pred), 2), n))

  r[n + 1,] <- 0
  N[,,n + 1] <- 0

  for (i in seq.int(n, 1))
  {
    ip1 <- i + 1

    if (!is.na(y[i]))
    {
      r[i,] <- tZ * v[i] / f[i] + t(L[,,i]) %*% r[ip1,]
      ahat[i,] <- a.pred[i,] + P.pred[,,i] %*% r[i,]

      if (i < convit || i > nmconvit)
      {
        N[,,i] <- crossprod(Z) / f[i] + crossprod(L[,,i], N[,,ip1]) %*% L[,,i]
        varahat[,,i] <- P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% P.pred[,,i]
      } else {
        N[,,i] <- N[,,ip1]
        varahat[,,i] <- varahat[,,ip1]
      }

    } else { # if (is.na(y[i]))
      # keep constant values, values from previous iteration
      # instead of NA since they are used in the next iteration
      r[i,] <- r[ip1,] #t(L[,,i]) %*% r[ip1,] # 'L' is NA
      ahat[i,] <- a.pred[i,] + P.pred[,,i] %*% r[i,]
      N[,,i] <- N[,,ip1]
      varahat[,,i] <- P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% P.pred[,,i]
    }

    ip2 <- i
    if (ip2 + 1 <= n && is.na(y[ip2]))
    { # iteration 'i+2' no longer used in this loop so it can be set to NA
      ip3 <- ip2 + 1 # the first row is removed below in 'r' and 'N'
      r[ip3,] <- NA
      ahat[ip2,] <- NA
      N[,,ip3] <- NA
      varahat[,,ip2] <- NA
    }
  }

  r <- matrix(r[-1,], nrow = nrow(r) - 1, ncol = ncol(r))
  N <- array(N[,,-1], dim = dim(N) - c(0,0,1))

  if (is.ts(y))
  {
    r <- ts(r)
    ahat <- ts(ahat)
    tsp(r) <- tsp(ahat) <- tsp(y)
  }    

  list(ahat = ahat, varahat = varahat, r = r, N = N)
}

KS.deriv <- function(y, ss, kf)
{
  n <- length(y)
  rp1 <- ncol(ss$V) + 1
  if (is.null(kf$convit)) {
    convit <- n + 1
    nmconvit <- n # no effect inthe if statement below
  } else {
    convit <- kf$convit
    nmconvit <- n - convit
  }

  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  tZ <- t(Z)
  tmT <- t(mT)

  v <- kf$v
  a.pred <- kf$a.pred
  P.pred <- kf$P.pred
  f <- kf$f
  L <- kf$L

  r <- matrix(nrow = n + 1, ncol = ncol(a.pred))
  N <- array(dim = c(ncol(a.pred), ncol(a.pred), n + 1))
  ahat <- matrix(nrow = n, ncol = ncol(a.pred))
  varahat <- array(NA, dim = c(rep(ncol(a.pred), 2), n))

  r[n + 1,] <- 0
  N[,,n + 1] <- 0

  da.pred <- kf$da.pred
  dP.pred <- kf$dP.pred
  df <- kf$df
  dL <- kf$dL
  dvof <- kf$dvof
  dr <- array(dim = c(n + 1, ncol(Z), rp1))
  dimnames(dr)[[3]] <- paste("var", seq(rp1), sep = "")
  dr[n + 1,,] <- 0
  dahat <- array(dim = c(n, ncol(Z), rp1))
  dN <- array(dim = c(n + 1, ncol(Z), ncol(Z), rp1))
  dimnames(dN)[[4]] <- dimnames(dr)[[3]]
  dN[n + 1,,,] <- 0
  dvarahat <- array(dim = c(n, ncol(Z), ncol(Z), rp1))

  fsq <- f^2
  ZtZ <- crossprod(Z)

  for (i in seq.int(n, 1))
  {
    ip1 <- i + 1

    if (!is.na(y[i]))
    {
      r[i,] <- tZ * v[i] / f[i] + t(L[,,i]) %*% r[ip1,]
      ahat[i,] <- a.pred[i,] + P.pred[,,i] %*% r[i,]

      if (i < convit || i > nmconvit)
      {
        N[,,i] <- ZtZ / f[i] + crossprod(L[,,i], N[,,ip1]) %*% L[,,i]
        varahat[,,i] <- P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% P.pred[,,i]
      } else {
        N[,,i] <- N[,,ip1]
        varahat[,,i] <- varahat[,,ip1]
      }
    } else { # if (is.na(y[i]))
      # keep constant values, values from previous iteration
      # instead of NA since they are used in the next iteration
      r[i,] <- r[ip1,] #t(L[,,i]) %*% r[ip1,] # 'L' is NA
      ahat[i,] <- a.pred[i,] + P.pred[,,i] %*% r[i,]
      N[,,i] <- N[,,ip1]
      varahat[,,i] <- P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% P.pred[,,i]
    }

    # derivative terms

    if (!is.na(y[i]))    
    {
      if (i < convit || i > nmconvit)
      {
        for (j in seq_len(rp1))
        {
          dr[i,,j] <- t(Z) * dvof[i,j] + t(dL[,,j,i]) %*% r[ip1,] +
            t(L[,,i]) %*% dr[ip1,,j]
          dN[i,,,j] <- -ZtZ * (df[i,j]/fsq[i]) + t(dL[,,j,i]) %*% N[,,ip1] %*% L[,,i] +
            t(L[,,i]) %*% dN[ip1,,,j] %*% L[,,i] + crossprod(L[,,i], N[,,ip1]) %*% dL[,,j,i]

          dahat[i,,j] <- (da.pred[i,,j] + dP.pred[,,i,j] %*% r[i,] + P.pred[,,i] %*% dr[i,,j])
          dvarahat[i,,,j] <- dP.pred[,,i,j] - dP.pred[,,i,j] %*% N[,,i] %*% P.pred[,,i] -
            P.pred[,,i] %*% dN[i,,,j] %*% P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% dP.pred[,,i,j]
        }
      } else { # if converged
        for (j in seq_len(rp1))
        {
          dr[i,,j] <- dr[ip1,,j]
          dN[i,,,j] <- dN[ip1,,,j]
          dahat[i,,j] <- dahat[ip1,,j]
          dvarahat[i,,,j] <- dvarahat[ip1,,,j]
        }
      }
    } else { # if (is.na(y[i]))
      for (j in seq_len(rp1))
      {
        dr[i,,j] <- dr[ip1,,j]
        dN[i,,,j] <- dN[ip1,,,j]
        dahat[i,,j] <- dahat[ip1,,j]
        dvarahat[i,,,j] <- dvarahat[ip1,,,j]
      }
    }

    ip2 <- i + 2
    if (ip2 <= n && is.na(y[ip2]))
    { # iteration 'i+2' no longer used in this loop so it can be set to NA      
      ip3 <- ip2 + 1 # the first row is removed below in 'r', 'N', 'dr' and 'dN'
      r[ip3,] <- NA
      ahat[ip2,] <- NA
      N[,,ip3] <- NA
      varahat[,,ip2] <- NA
      dr[ip3,,] <- 0
      dN[ip3,,,] <- 0
      dahat[ip2,,] <- 0
      dvarahat[ip2,,,] <- 0
    }
  }

  r <- ts(matrix(r[-1,], nrow = nrow(r)-1, ncol = ncol(r)), 
    frequency = frequency(y), end = end(y))
  N <- array(N[,,-1], dim = dim(N) - c(0,0,1))

  # in the local level model the second dimension is lost
  if (ncol(Z) > 1)
  {
    dr <- dr[-1,,]
    dN <- dN[-1,,,]
  } else {
    dr <- array(dr[-1,,], dim = c(n, 1, rp1))
    dN <- array(dN[-1,,,], dim = c(n, 1, 1, rp1))
  }

  ahat <- ts(ahat)
  tsp(ahat) <- tsp(y)

  list(ahat = ahat, varahat = varahat, r = r, N = N,
    dahat = dahat, dvarahat = dvarahat, dr = dr, dN = dN)
}
