#' Print method for objects of class \code{Bolstad}
#' 
#' This function provides a print summary method for the output of
#' \code{bayes.lm}.
#' 
#' 
#' 
#' @param x an object of class \code{Bolstad}
#' @param digits number of digits to print 
#' @param \dots any other arguments that are to be passed to \code{print.default}
#' @details if x has both class \code{Bolstad} and \code{lm} then a print method 
#' similar to \code{print.lm} is called, otherwise \code{print.default} is called
#' @author James Curran
#' @seealso \code{\link{bayes.lm}}
#' @export

print.Bolstad = function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if(length(class(x)) == 2 && all(grepl("Bolstad|lm", class(x)))){
    getTermLabels = function(x) 
      attr(x$terms, "term.labels")
    
    cat("\nCall:", paste0(deparse(x$call)), sep = "\n", collapse= "\n")
    
    if (length(coef(x))) {
      cat("Coefficients:\n")
      c = coef(x)
      names(c) = c("(Intercept)", getTermLabels(x))
      print(format(c, digits = digits), print.gap = 2L, quote = FALSE)
    }
    cat("\n")
  }else{
    print.default(x, ...)
  }
}

