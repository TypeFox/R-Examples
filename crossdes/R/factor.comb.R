factor.comb <- function(p, n) {
  f <- matrix(0, p^n, n)
  for (i in 1:n) {
    f[, i] <- rep(0:(p - 1), p^(i - 1), each = p^(n - i))
  }
  f
} 
