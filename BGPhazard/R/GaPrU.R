GaPrU <-
function(alpha, beta, c.r, lambda.r) {
  uk <- 100
  K <- length(lambda.r)
  matfu <- matrix(0, nrow = (uk + 1), ncol = (K - 1))
  matfu[1, ] <- 1
  for(k in 1:(K - 1)) {
    if (c.r[k] != 0) {
      for(u in 0:uk) {
        logmatfu <- (u * (log(c.r[k]) + log(c.r[k] + beta[k + 1]) 
                          + log(lambda.r[k]) + log(lambda.r[k + 1])) 
                     - lgamma(u + 1) - lgamma(alpha[k + 1] + u))
        matfu[u + 1, k] <- exp(logmatfu)
      }
    }
  }
  matfu <- prop.table(matfu, margin = 2)
  return(matfu)
}
