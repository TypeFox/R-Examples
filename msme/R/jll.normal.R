jll_normal <- function(params, X, y) {
  p <- length(params)
  beta <- params[-p]
  sigma <- exp(params[p])
  linpred <- X %*% beta  
  sum(dnorm(y, mean = linpred, sd = sigma, log = TRUE))
}
