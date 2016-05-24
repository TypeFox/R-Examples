BePrC <-
function(alpha, beta, Pi.r, u.r, epsilon) {
  K <- length(Pi.r)
  ck <- 50
  matfc <- matrix(0, nrow = (ck + 1), ncol = (K - 1))
  for(k in 1:(K - 1)) {
    for(i in u.r[k]:(u.r[k] + ck)) {
      lgammas <- (lgamma(alpha[k + 1] + beta[k + 1] + i) 
                  - lgamma(beta[k + 1] + i - u.r[k]) - lgamma(i - u.r[k] + 1))
      logmatfc <- (lgammas + i * (log(epsilon) + log(1 - Pi.r[k]) 
                                  + log(1 - Pi.r[k + 1])))
      matfc[i - u.r[k] + 1, k] <- exp(logmatfc)
    }
  }
  matfc <- prop.table(matfc, margin = 2)
  return(matfc)
}
