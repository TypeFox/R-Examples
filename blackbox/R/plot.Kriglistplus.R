plot.OKriglistplus <- function(x, digits = 4, which = 1, ...) {
  zut <- x[c("y", "fitted.values")]
  class(zut) <- c("OKrig") ## not a true OKrig object, but good enough for call to plot.OKrig
  plot(zut, ...)
}
