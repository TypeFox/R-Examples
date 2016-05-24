
predict.stsmSS <- function(object, y, n.ahead = 12L, ...)
{
  kf <- KF(y, object)
  
  n <- length(y)
  Z <- object$Z
  mT <- object$T
  H <- object$H
  Q <- object$Q
  tmT <- t(mT)

  m <- ncol(Z)
  a <- matrix(nrow = n.ahead + 1, ncol = m)
  P <- array(NA, dim = c(m, m, n.ahead + 1))
  pred <- var <- rep(0, n.ahead)

  a[1,] <- if (m == 1) kf$a.upd[n] else kf$a.upd[n,]
  P[,,1] <- if (m == 1) kf$P.upd[n] else kf$P.upd[,,n]

  for (i in seq_len(n.ahead))
  {
    ip1 <- i + 1
    a[ip1,] <- mT %*% a[i,]
    P[,,ip1] <- mT %*% P[,,i] %*% tmT + Q

    pred[i] <- Z %*% a[ip1,]
    var[i] <- H + sum(diag(crossprod(Z) %*% P[,,ip1]))
  }

  s <- frequency(y)
  ytsp <- tsp(y)
  t0 <- ytsp[2L] + 1 / s

  a <- ts(a[-1,], start = t0, frequency = s)
  Pres <- ts(matrix(nrow = n.ahead, ncol = m), start = t0, frequency =  s)
  if (m == 1)
  {
   Pres <- ts(P[-1], start = t0, frequency =  s)
  } else {
    P <- P[,,-1]
    for (i in seq_len(n.ahead))
      Pres[i,] <- diag(P[,,i])
  }

  pred <- ts(pred, start = t0, frequency =  s)
  var <- ts(var, start = t0, frequency =  s)

  list(pred = pred, se = sqrt(var), a = a, P = Pres)
}
