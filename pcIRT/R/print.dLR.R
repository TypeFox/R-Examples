#' @rdname dLR
#'
#' @export
#'
#' @method print dLR
#'
#' @param x object of class \code{dLR}

print.dLR <-
function(x, ...){
    
    cat("\n Call: \t Likelihood Ratio Test MPRM: dimension reduction \n\n")
    
    cat("emp Chi2: \t", x$emp_Chi2, "\n")    
    cat("df: \t\t\t\t\t\t\t", x$df, "\n")
    cat("p-value: \t\t", deparse(round(x$pval,3)), "\n\n")
    
    
  }
