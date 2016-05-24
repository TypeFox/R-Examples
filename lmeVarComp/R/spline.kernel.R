spline.kernel <-
function(x, normalized = FALSE)
{
  n <- nrow(x)
  k <- seq_len(n)
  i <- rep(k, times = n)
  j <- rep(k, each = n)
  xi <- x[i, , drop = FALSE]
  xj <- x[j, , drop = FALSE]
  xm <- pmin(xi, xj)
  u <- 1 + xi * xj + xi * xj * xm / 2 - xm ^ 3 / 6
  if (normalized) {
    v1 <- 1 + xi * xi + xi * xi * xi / 3
    v2 <- 1 + xj * xj + xj * xj * xj / 3
    u <- u / sqrt(v1 * v2)
  }
  kernel <- matrix(exp(rowSums(log(u))), n, n)
  kernel
}
