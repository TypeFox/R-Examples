plot.rtif <- function(x, theta = x$theta,
  xlab = expression(theta), ylab = "Information",
  main = "TIF", lwd = 2, ...) {

  plot.default(theta, x$pq, xlab = xlab, ylab = ylab,
    main = main, type = "l", lwd = lwd, ...)
}
