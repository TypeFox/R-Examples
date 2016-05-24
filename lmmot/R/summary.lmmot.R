##' Print details about lmmot object.
##'
##' @title Summary if lmmot object
##'
##' @param object lmmot object to print.
##' @param digits number of decimal digits to print.
##' @param ... further arguments passed to or from other methods.
##'
##' @export
##' @seealso \link[stats]{lm} \link[lmmot]{lmmot}
##' @author Marvin Wright

summary.lmmot <- function (object, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  
  cat("\nCensoring:\n")
  print(object$censoring)
  cat("\n")
  
  if (length(coef(object))) {
    cat("Coefficients:\n")
    
    res <- cbind("Estimate"=object$estimate,
                 "Std. error"=object$stdEr,
                 "t value"=object$tval, 
                 "Pr(> t)"=object$pval)
    
    print.default(format(res, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  
  invisible(object)
}