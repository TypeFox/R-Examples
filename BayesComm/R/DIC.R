DIC <- function (BC) {  # make this a 
  Z <- BC$trace$z
  Y <- BC$call$Y
  # calculate D_hat
  zhat <- apply(Z, c(2, 3), mean)  # posterior mean of z
  p <- pnorm(zhat, log.p = TRUE)
  ind <- which(!Y)
  p[ind] <- log(-expm1(p[ind]))
  Dhat <- -2 * sum(p)
  # calculate D_bar
  p <- pnorm(Z, log.p = TRUE)
  for(j in 1:ncol(Y)) {
    ind <- which(!Y[, j])
    p[, ind, j] <- log(-expm1(p[, ind, j]))
  }
  Dbar <- mean(-2 * apply(p, 1, sum))
  # return DIC
  2 * sum(Dbar) - Dhat
}