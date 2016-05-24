print.lm.beta <- function(x, standardized=TRUE, ...) {
  if(standardized) {
    cat("\nCall:\n")
    print(x$call,...)
    cat("\nStandardized Coefficients::\n")
    print(x$standardized.coefficients,...)
    cat("\n")
  } else {
    x.tmp <- x
    attr(x.tmp,"class") <- "lm"
    print(x.tmp,...)
  }
  invisible(x)
}