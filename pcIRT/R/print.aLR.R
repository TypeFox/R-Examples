#' @rdname lrt
#'
#' @export
#'
#' @method print aLR
#' @param x Object of class aLR
#' @param \ldots further arguments



print.aLR <-
function(x,...){
    
    cat("\n Call: \t Likelihood Ratio Test MPRM \n\n")
    
    cat("emp Chi2: \t", x$emp_Chi2, "\n")    
    cat("df: \t\t\t\t\t\t\t", x$df, "\n")
    cat("p-value: \t\t", deparse(round(x$pval,3)), "\n\n")    
  }
