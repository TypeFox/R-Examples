print.tpfit <-
function(x, ...) {
  cat("Estimated transition rates:\n\n")
  print(unclass(x$coefficients), ...)
  cat("\n")
  if (!is.null(x$prop)) {
    cat("Estimated proportions:\n\n")
    proportions <- as.vector(x$prop)
    names(proportions) <- names(x$prop)
    print(proportions, ...)
    cat("\n")
  }
  invisible(x)
}
