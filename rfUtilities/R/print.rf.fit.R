print.rf.fit <- function(x, ...) {
  cat(x$message, "(Ideal overfit ratio >= 10)", "\n")
  cat("", "\n")
  cat("Fit statistics", "\n")
  print(x$fit)
}
