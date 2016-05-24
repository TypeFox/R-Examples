#' @rdname dLR
#'
#' @export
#'
#' @method summary dLR
#'
#' @param object object of class \code{dLR}
#' @param \dots \dots{}

summary.dLR <-
function(object, ...){
    
    cat("\n Call: \t Likelihood Ratio Test MPRM: dimension reduction \n\n")
    
    cat("emp Chi2: \t", object$emp_Chi2, "\n")    
    cat("df: \t\t\t\t\t\t\t", object$df, "\n")
    cat("p-value: \t\t", deparse(round(object$pval,3)), "\n\n")
    
    
  }
