print.minmax <- function (x, digits=max(3L, getOption("digits") - 3L), ...) {
  cat("Minima: \n")
  print(x$mins, digits=digits)
  cat("\nMaxima: \n")
  print(x$maxs, digits=digits)
}
