plot.forensicbatwing <-
function(x, ...) {
  plot(x$result[, 1], x$result[, 2], type = "l", main = "Traceplot", xlab = "Iteration", ylab = "Estimate of match probability")
  invisible(NULL)
}
