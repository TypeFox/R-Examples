UpdPi <-
function(alpha, beta, c.r, u.r, n, m) {
  K <- length(alpha)
  Pi.r <- rep(0, K)
  Pi.r[1] <- rbeta(1, shape1 = alpha[1] + u.r[1] + n[1], 
                   shape2 = beta[1] + c.r[1] - u.r[1] + m[1])
  for(k in 2:(K - 1)) {
    a <- alpha[k] + u.r[k - 1] + u.r[k] + n[k]
    b <- beta[k] + c.r[k - 1] - u.r[k - 1] + c.r[k] - u.r[k] + m[k]
    if (a > 0 && b > 0) {
      Pi.r[k] <- rbeta(1, shape1 = a, shape2 = b)
    }
  }
  Pi.r[K] <- rbeta(1, shape1 = alpha[K] + u.r[K - 1] + n[K], 
                   shape2 = beta[K] + c.r[K - 1] - u.r[K - 1] + m[K])
  for (k in 1:K) {
    Pi.r[k] <- min(max(0.001, Pi.r[k]), 0.999)
  }
  return(Pi.r)
}
