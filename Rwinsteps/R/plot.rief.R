plot.rief <- function(x, theta = x$theta, add = FALSE,
  xlab = expression(theta), ylab = "Standard Error", main = "IEF",
  lwd = 2, col = "r", ...) {

  ni <- ncol(x$p)

  if(!add)
    plot.default(c(min(theta), max(theta)), c(0, max(x$e)),
      xlab = xlab, ylab = ylab, main = main, type = "n", ...)

  lines(x, theta, lwd = lwd, col = col)
}
