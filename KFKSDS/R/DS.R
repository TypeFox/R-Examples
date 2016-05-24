
DS <- function(y, ss, kf, ks)
{
  n <- length(y)
  if (is.null(kf$convit)) {
    convit <- n + 1
    nmconvit <- n
  } else {
    convit <- kf$convit
    nmconvit <- n - convit
  }

  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  R <- ss$R
  V <- ss$V
  tZ <- t(Z)
  tmT <- t(mT)
  VtR <- tcrossprod(V, R)

  v <- kf$v
  f <- kf$f
  K <- kf$K

  r <- ks$r
  N <- ks$N

  epshat <- rep(NA, n)
  vareps <- rep(NA, n)
  etahat <- matrix(nrow = n, ncol = ncol(ss$V))
  vareta <- array(dim = c(n, ncol(ss$V), ncol(ss$V)))

  for (i in seq(n, 1))
  {
    cKi <- cbind(K[i,])
    if (!is.na(y[i]))
    {
      epshat[i] <- H * (v[i] / f[i] - crossprod(cKi, r[i,]))
      etahat[i,] <- VtR %*% cbind(r[i,])

      if (i < convit || i > nmconvit)
      {
        vareps[i] <- H - H^2 * (1/f[i] + crossprod(cKi, N[,,i]) %*% cKi)
        vareta[i,,] <- V - VtR %*% N[,,i] %*% R %*% V
      } else {
        vareps[i] <- vareps[i+1]
        vareta[i,,] <- vareta[i+1,,]
      }
    } else { # if (is.na(y[i]))
      #epshat[i] <- NA
      #vareta[i,,] <- NA
    }
  }

  if (is.ts(y))
  {
    epshat <- ts(epshat)
    etahat <- ts(etahat)
    tsp(epshat) <- tsp(etahat) <- tsp(y)
  }

  list(epshat = epshat, vareps = vareps, 
    etahat = etahat, vareta = vareta)
}

DS.deriv <- function(ss, ksd)
{
  n <- nrow(ksd$ahat)
  rp1 <- ncol(ss$V) + 1
  dahat <- ksd$dahat
  dvarahat <- ksd$dvarahat
  r <- ksd$r
  N <- ksd$N
  dr <- ksd$dr
  dN <- ksd$dN

  depshat <- dvareps <- rep(NA, rp1)
  ref <- which(ss$Z == 1)
  for (j in seq_along(depshat))
  {
    depshat[j] <- -sum(dahat[,ref,j])
    dvareps[j] <- sum(dvarahat[,ref,ref,j])
  }

  ref <- which(diag(ss$R) == 1)
  detahat <- matrix(0, nrow = length(ref), ncol = rp1)
  for (id in seq_along(ref))
  {
    detahat[id,] <- ss$V[id,id] * colSums(dr[,ref[id],])
    detahat[id,id+1] <- detahat[id,id+1] + sum(r[,ref[id]])
  }

  dvareta <- matrix(0, nrow = length(ref), ncol = rp1)
  for (id in seq_along(ref))
  {
    V <- ss$V[id,id] #diag(ss$V)[id]
    rid <- ref[id]
    dvareta[id,] <- -V^2 * colSums(dN[,rid,rid,])
    dvareta[id,id+1] <- dvareta[id,id+1] + n - 2 * V * sum(N[rid,rid,])
  }

  list(depshat = depshat, dvareps = dvareps, 
    detahat = detahat, dvareta = dvareta)
}
