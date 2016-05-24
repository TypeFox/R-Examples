AICc <- function(obj) {
  N <- length(obj$time)
  P <- obj$Pn
  SSE <- obj$val
  AIC <- N * log(SSE / N) + 2 * P
  AIC + (2 * P * (P + 1)) / (N - P - 1)
}
