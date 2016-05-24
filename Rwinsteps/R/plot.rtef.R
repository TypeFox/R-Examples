plot.rtef <- function(x, theta = x$theta,
  xlab = expression(theta), ylab = "Standard Error",
  main = "TEF", lwd = 2, ...) {

  plot.default(theta, x$e, xlab = xlab, ylab = ylab,
    main = main, type = "l", lwd = lwd, ...)
}
