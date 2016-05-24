print.oneway.rfit <- function (x, digits = max(5, .Options$digits - 2), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nOverall Test of All Locations Equal \n")
  print(drop.test(x$fit))
  cat("\n")
  print(x$pp)
}
