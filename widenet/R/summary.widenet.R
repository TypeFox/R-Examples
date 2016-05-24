## This file has summary, print.summary, and print methods for widenet objects

summary.widenet <- function(object, ...) {

  summary.result <- object[c("call", "min.cvm.mat", "total.time")]

  class(summary.result) <- "summary.widenet"
  
  return(summary.result)
}

print.summary.widenet <- function(x, digits = 3, ...) {

  cat("The call was:\n\n", deparse(x$call), "\n\n")

  cat("The minimum cross-validated risk by order and alpha value:\n\n")

  print(x$min.cvm.mat, ...)

  cat("\nThe total time taken to fit this model was:\n\n")

  print(x$total.time, digits, ...)

  cat("\nseconds.\n")

  invisible(x)
}

print.widenet <- function(x, digits = 3, ...) {

  print(summary(x), digits, ...)

  invisible(x)
}
