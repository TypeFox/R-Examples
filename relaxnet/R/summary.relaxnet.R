## This file has summary, print.summary, and print methods for relaxnet objects

summary.relaxnet <- function(object, ...) {

  summary.result <- object[c("call", "relax.num.vars", "total.time")]

  class(summary.result) <- "summary.relaxnet"
  
  return(summary.result)
}

print.summary.relaxnet <- function(x, digits = 3, ...) {

  cat("The call was:\n\n", deparse(x$call), "\n\n")

  cat("Relaxed models were fit containing the following number of",
      "variables:\n\n")

  print(x$relax.num.vars, ...)

  cat("\nThe total time taken to fit this model was:\n\n")

  print(x$total.time, digits, ...)

  cat("\nseconds.\n")

  invisible(x)
}

print.relaxnet <- function(x, digits = 3, ...) {

  print(summary(x), digits, ...)

  invisible(x)
}
