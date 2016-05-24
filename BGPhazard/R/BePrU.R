BePrU <-
function(alpha, beta, c.r, Pi.r) {
  K <- length(Pi.r)
  matfu <- matrix(0, nrow = (max(c.r) + 1), ncol = (K - 1))
  matfu[1, ] <- 1
  for(k in 1:(K - 1)) {
    if (c.r[k] != 0) {
      for(u in 0:c.r[k]) {
        lphi <- (log(Pi.r[k]) + log(Pi.r[k + 1]) 
                 - log(1 - Pi.r[k]) - log(1 - Pi.r[k + 1]))
        logmatfu <- (u * lphi - lgamma(u + 1) - lgamma(c.r[k] - u + 1) 
                     - lgamma(alpha[k + 1] + u) 
                     - lgamma(beta[k + 1] + c.r[k] - u))
        matfu[u + 1,  k] <- exp(logmatfu)
      }
    }
  }
  matfu <- prop.table(matfu, margin = 2)
  return(matfu)
}
