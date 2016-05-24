
KSDS <- function(y, ss, kf)
{
  n <- length(y)

  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  tZ <- t(Z)
  tmT <- t(mT)
  R <- ss$R
  V <- ss$V

  v <- kf$v
  a.pred <- kf$a.pred
  P.pred <- kf$P.pred
  f <- kf$f
  K <- kf$K
  L <- kf$L

  r <- matrix(nrow = n + 1, ncol = ncol(a.pred))
  N <- array(dim = c(ncol(a.pred), ncol(a.pred), n + 1))
  ahat <- matrix(nrow = n, ncol = ncol(a.pred))
  varahat <- array(NA, dim = c(rep(ncol(a.pred), 2), n))

  r[n + 1,] <- 0
  N[,,n + 1] <- 0

  epshat <- rep(NA, n)
  vareps <- matrix(nrow = n, ncol = 1)
  etahat <- matrix(nrow = n, ncol = ncol(ss$V))
  vareta <- array(dim = c(n, ncol(ss$V), ncol(ss$V)))

  for (i in seq(n, 1))
  {
    ip1 <- i + 1
    r[i,] <- tZ * v[i] / f[i] + t(L[,,i]) %*% r[ip1,]
    N[,,i] <- crossprod(Z) / f[i] + crossprod(L[,,i], N[,,ip1]) %*% L[,,i]

    ahat[i,] <- a.pred[i,] + P.pred[,,i] %*% r[i,]
    varahat[,,i] <- P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% P.pred[,,i]

    # disturbance smoother

    epshat[i] <- H * (v[i] / f[i] - crossprod(cbind(K[i,]), r[ip1,]))
    vareps[i,] <- H - H^2 * (1/f[i] + crossprod(cbind(K[i,]), N[,,ip1]) %*% cbind(K[i,]))

    tmp <- tcrossprod(V, R)
    etahat[i,] <- tmp %*% cbind(r[ip1,])
    vareta[i,,] <- V - tmp %*% N[,,ip1] %*% R %*% V
  }

  r <- ts(matrix(r[-1,], nrow = nrow(r)-1, ncol = ncol(r)), 
    frequency = frequency(y), end = end(y))
  N <- array(N[,,-1], dim =  dim(N) - c(0,0,1))

  ahat <- ts(ahat)
  tsp(ahat) <- tsp(y)

  epshat <- ts(epshat)
  etahat <- ts(etahat)
  tsp(epshat) <- tsp(etahat) <- tsp(y)

  list(epshat = epshat, vareps = vareps, 
    etahat = etahat, vareta = vareta,
    r = r, N = N, ahat = ahat, varahat = varahat)
}

KSDS.deriv <- function(y, ss, kf)
{
  n <- length(y)

  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  tZ <- t(Z)
  tmT <- t(mT)
  R <- ss$R
  V <- ss$V
  tmT <- t(mT)

  v <- kf$v
  a.pred <- kf$a.pred
  P.pred <- kf$P.pred
  f <- kf$f
  L <- kf$L
  K <- kf$K

  r <- matrix(nrow = n + 1, ncol = ncol(a.pred))
  N <- array(dim = c(ncol(a.pred), ncol(a.pred), n + 1))
  ahat <- matrix(nrow = n, ncol = ncol(a.pred))
  varahat <- array(NA, dim = c(rep(ncol(a.pred), 2), n))

  r[n + 1,] <- 0
  N[,,n + 1] <- 0

nv <- ncol(ss$V) + 1

da.pred <- kf$da.pred
dP.pred <- kf$dP.pred
df <- kf$df
dL <- kf$dL
dvof <- kf$dvof
dr <- array(dim = c(n + 1, ncol(Z), nv))
dimnames(dr)[[3]] <- paste("var", seq(nv), sep = "")
dr[n + 1,,] <- 0
dahat <- array(dim = c(n, ncol(Z), nv))
dN <- array(dim = c(n + 1, ncol(Z), ncol(Z), nv))
dimnames(dN)[[4]] <- dimnames(dr)[[3]]
dN[n + 1,,,] <- 0
dvarahat <- array(dim = c(n, ncol(Z), ncol(Z), nv))

fsq <- f^2

ZtZ <- crossprod(Z)

  epshat <- rep(NA, n)
  vareps <- rep(NA, n)
  etahat <- matrix(nrow = n, ncol = ncol(ss$V))
  vareta <- array(dim = c(n, ncol(ss$V), ncol(ss$V)))

  for (i in seq(n, 1))
  {
    ip1 <- i + 1
    r[i,] <- tZ * v[i] / f[i] + t(L[,,i]) %*% r[ip1,]
    N[,,i] <- ZtZ / f[i] + crossprod(L[,,i], N[,,ip1]) %*% L[,,i]

    ahat[i,] <- a.pred[i,] + P.pred[,,i] %*% r[i,]
    varahat[,,i] <- P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% P.pred[,,i]

for (j in seq(nv))
{
dr[i,,j] <- t(Z) * dvof[i,j] + t(dL[,,j,i]) %*% r[ip1,] +
  t(L[,,i]) %*% dr[ip1,,j]
dN[i,,,j] <- -ZtZ * (df[i,j]/fsq[i]) + t(dL[,,j,i]) %*% N[,,ip1] %*% L[,,i] +
  t(L[,,i]) %*% dN[ip1,,,j] %*% L[,,i] + crossprod(L[,,i], N[,,ip1]) %*% dL[,,j,i]

dahat[i,,j] <- (da.pred[i,,j] + dP.pred[,,i,j] %*% r[i,] + P.pred[,,i] %*% dr[i,,j])
dvarahat[i,,,j] <- dP.pred[,,i,j] - dP.pred[,,i,j] %*% N[,,i] %*% P.pred[,,i] -
  P.pred[,,i] %*% dN[i,,,j] %*% P.pred[,,i] - P.pred[,,i] %*% N[,,i] %*% dP.pred[,,i,j]
}

    # Disturbance smoother

    epshat[i] <- H * (v[i] / f[i] - crossprod(cbind(K[i,]), r[ip1,]))
    vareps[i] <- H - H^2 * (1/f[i] + crossprod(cbind(K[i,]), N[,,ip1]) %*% cbind(K[i,]))

    tmp <- tcrossprod(V, R)
    etahat[i,] <- tmp %*% cbind(r[ip1,])
    vareta[i,,] <- V - tmp %*% N[,,ip1] %*% R %*% V
  }

  r <- ts(matrix(r[-1,], nrow = nrow(r)-1, ncol = ncol(r)), 
    frequency = frequency(y), end = end(y))
  N <- array(N[,,-1], dim =  dim(N) - c(0,0,1))

if (ncol(Z) > 1)
{
dr <- dr[-1,,]
dN <- dN[-1,,,]
} else
{
dr <- array(dr[-1,,], dim = c(n, 1, nv))
dN <- array(dN[-1,,,], dim = c(n, 1, 1, nv))
}

  ahat <- ts(ahat)
  tsp(ahat) <- tsp(y)

  epshat <- ts(epshat)
  etahat <- ts(etahat)
  tsp(epshat) <- tsp(etahat) <- tsp(y)

  list(r = r, N = N, ahat = ahat, varahat = varahat, 
    dr = dr, #dr = colSums(dr), 
    dN = dN, #dN = colSums(dN)
    dahat = dahat, dvarahat = dvarahat,
    epshat = epshat, vareps = vareps, 
    etahat = etahat, vareta = vareta)
}
