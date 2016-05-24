# S^2 estimation
s2 <- function(y, w = rep(1, length(y))) {
  N <- sum(w)
  n <- length(y)
  s2 <- (N - 1) / N * n / (n - 1) *
    (sum(y ^ 2 * w) - sum(y * w) ^ 2 / N) / (N - 1)
  return(s2)
}
