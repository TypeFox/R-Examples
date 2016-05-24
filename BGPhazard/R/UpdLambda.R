UpdLambda <-
function(alpha, beta, c.r, u.r, n, m) {
  K <- length(alpha)
  lambda.r <- rep(0, K)
  while(lambda.r[1] == 0) {
    lambda.r[1] <- rgamma(1, shape = alpha[1] + u.r[1] + n[1], 
                          scale = 1 / (beta[1] + c.r[1] + m[1]))
  }
  for(k in 2:(K - 1)) {
    while(lambda.r[k] == 0) {
      lambda.r[k] <- rgamma(1, shape = alpha[k] + u.r[k - 1] + u.r[k] + n[k], 
                            scale = 1 / (beta[k] + c.r[k - 1] + c.r[k] + m[k]))
    }
  }
  while(lambda.r[K] == 0) {
    lambda.r[K] <- rgamma(1, shape = alpha[K] + u.r[K - 1] + n[K], 
                          scale = 1 / (beta[K] + c.r[K - 1] + m[K]))
  }
  return(lambda.r)
}
