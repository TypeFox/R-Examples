GaUpdU <-
function(alpha, beta, c.r, lambda.r) {
  f.umat <- GaPrU(alpha, beta, c.r, lambda.r)
  K <- length(lambda.r)
  u.r <- rep(0, K - 1)
  for(i in 1:(K - 1)) {
    u.r[i] <- sample(x = 0:100, size = 1, prob = f.umat[, i])
  }
  return(u.r)
}
