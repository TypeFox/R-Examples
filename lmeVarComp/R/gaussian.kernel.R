gaussian.kernel <-
function(x, rho = 10)
{
  n <- nrow(x)
  k <- seq_len(n)
  i <- rep(k, times = n)
  j <- rep(k, each = n)
  u <- x[i, , drop = FALSE] - x[j, , drop = FALSE]
  kernel <- matrix(exp(-rho * rowSums(u ^ 2)), n, n)
  kernel
}
