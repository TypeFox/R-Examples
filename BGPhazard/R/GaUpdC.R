GaUpdC <-
function(alpha, beta, c.r, lambda.r, u.r, epsilon) {
  K <- length(lambda.r)
  c.str <- rep(0, K - 1)
  for(k in 1:(K - 1)) {
    c.str[k] <- rgamma(1, shape = u.r[k] + 1, 
                       scale = 1 / (lambda.r[k + 1] + lambda.r[k] + epsilon))
  }
  for(k in 1:(K - 1)) {
    lw.c.str <- (alpha[k + 1] + u.r[k]) * log(beta[k + 1] + c.str[k])
    lw.c.r <- (alpha[k + 1] + u.r[k]) * log(beta[k + 1] + c.r[k])
    ratio <- lw.c.str - lw.c.r
    p <- min(exp(ratio), 1)
    if (runif(1) <= p) {
      c.r[k] <- c.str[k]
    }
  }
  return(c.r)
}
