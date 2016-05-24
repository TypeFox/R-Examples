plot.rtrf <- function(x, theta = x$theta,
  xlab = expression(theta), ylab = "Total", main = "TRF",
  lwd = 2, ...) {

  plot.default(theta, x$p, xlab = xlab, ylab = ylab,
    main = main, type = "l", lwd = lwd, ...)
}
