##' Print lmmot object.
##'
##' @title Print lmmot object
##'
##' @param x lmmot object to print.
##' @param digits number of decimal digits to print.
##' @param ... further arguments passed to or from other methods.
##'
##' @export
##' @seealso \link[stats]{lm} \link[lmmot]{lmmot}
##' @author Marvin Wright

print.lmmot <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  
  cat("\nCensoring:\n")
  print(x$censoring)
  cat("\n")
  
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}