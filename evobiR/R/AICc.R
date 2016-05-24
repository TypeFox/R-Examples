AICc <- function(loglik, K, N){
  return(-2 * loglik + 2 * K * (N / ( N - K - 1)))
}
