plot.riif <- function(x, theta = x$theta, add = FALSE,
  xlab = expression(theta), ylab = "Information", main = "IIF",
  lwd = 2, col = "r", ...) {

  ni <- ncol(x$pq)

  if(!add)
    plot.default(c(min(theta), max(theta)), c(0, max(x$pq)),
      xlab = xlab, ylab = ylab, main = main, type = "n", ...)

  lines(x, theta, lwd = lwd, col = col)
}
