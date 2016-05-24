#' @rdname lrt
#'
#' @export
#'
#' @method summary aLR
#' 

summary.aLR <-
function(object,...){
    
    cat("\n Call: \t Likelihood Ratio Test MPRM \n\n")
    
    cat("emp Chi2: \t", object$emp_Chi2, "\n")    
    cat("df: \t\t\t\t\t\t\t", object$df, "\n")
    cat("p-value: \t\t", deparse(round(object$pval,3)), "\n\n")
    
    cat("--------------------------------------------------------------------\n")
    cat("Parameter estimates: \n")
    print(object$itempar)
    cat("SE estimates: \n")
    print(object$item_se)
    
  }
