BeUpdU <-
function(alpha,  beta,  c.r,  Pi.r) {
  f.umat <- BePrU(alpha,  beta,  c.r,  Pi.r)
  K <- length(Pi.r)
  u.r <- rep(0, K - 1)
  for(i in 1:(K - 1)) {
    u.r[i] <- sample(x = 0:max(c.r), size = 1, prob = f.umat[, i])
  }
  return(u.r)
}
