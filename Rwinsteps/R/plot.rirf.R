plot.rirf <- function(x, theta = x$theta, add = FALSE,
  xlab = expression(theta), ylab = "P(X)", main = "IRF",
  lwd = 2, col = "r", ...) {

  ni <- ncol(x$p)

  if(!add)
    plot.default(c(min(theta), max(theta)), c(0, 1),
      xlab = xlab, ylab = ylab, main = main, type = "n", ...)

  lines(x, theta, lwd = lwd, col = col)
}
