sample_e <-
function(e, trunc, iR) {
  e2 <- e
  n <- dim(e)[1]
  nsp <- dim(iR)[1]
  coy2 <- 1 / iR[col(iR) == row(iR)]
  A <- iR
  A[col(A) == row(A)] <- 0
  for (i in 1:nsp) {
    mn <- -coy2[i] * (e2 %*% A[i, ])
    std <- sqrt(coy2[i])
    n <- length(mn)
    e2[, i] <- rtnorm(n, mn, rep(std, n), trunc[, i, 1], trunc[, i, 2])
  }
  e2
}