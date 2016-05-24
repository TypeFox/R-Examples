BeUpdC <-
function(alpha, beta, Pi.r, u.r, epsilon) {
  f.cmat <- BePrC(alpha, beta, Pi.r, u.r, epsilon)
  K <- length(Pi.r)
  c.r <- rep(0, K - 1)
  ck <- 50
  for(k in 1:(K - 1)) {
    c.r[k] <- sample(x = u.r[k]:(u.r[k] + ck), size = 1, prob = f.cmat[, k])
  }
  return(c.r)
}
