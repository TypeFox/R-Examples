print.eiReg <- function(x, digits = max(2, getOption("digits") - 4), ...) { 
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("Estimated internal cells:\n")
  print.default(format(x$coef, digits = digits), 
                print.gap = 2, quote = FALSE, ...)
}
