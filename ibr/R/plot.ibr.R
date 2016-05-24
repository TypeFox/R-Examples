plot.ibr <- function(x,...) {
  r <- residuals(x)
  yh <- predict(x)
  plot(yh, r, xlab = "Fitted", ylab = "Residuals",...)
  invisible()
}
