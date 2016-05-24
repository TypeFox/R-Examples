print.endogMNP <- function (x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    param <- apply(x$param, 2, mean)
    if (length(param)) {
      cat("Parameter estimates (posterior means):\n")
      print.default(format(param, digits = digits), print.gap = 2,
                    quote = FALSE)
    }
    else cat("No parameter estimates\n")
    cat("\n")
    invisible(x)
  }
