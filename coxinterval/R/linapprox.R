### linear approximation
linapprox <- function(xyin, xout)
{
  x <- sort(unique(c(0, xyin[, 1])))
  y <- as.vector(tapply(c(0, xyin[, 2]), c(0, xyin[, 1]), max))
  n <- length(x)
  xmax <- max(xout)
  if (max(x) < xmax) {
    x <- c(x, xmax)
    y <- c(y, y[n-1] + (xmax - x[n-1])/(x[n] - x[n-1]) * (y[n] - y[n-1]))
  }
  approx(x, y, xout)$y
}
