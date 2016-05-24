print.aspect <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nValue target function:",round(x$loss,3))
  cat("\n\nCorrelation matrix of the transformed data:\n")
  print(round(x$cormat, 3))
  cat("\n")
  invisible(x)
}
