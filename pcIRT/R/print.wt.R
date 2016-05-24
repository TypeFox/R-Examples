#' @rdname wt
#'
#' @export
#'
#' @method print wt
#'
#' @param x object of class \code{wt}


print.wt <-
function(x, ...){
  
  cat("\n Weight Test unidimensional polytomous model, \n\n")
    
  cat("Score parameters: \n")
  print(x$unconstrScoreParameter)
  
  if(!is.null(x$emp_Chi2)){
    cat("-----------------------------------------------------\n")
        
    cat("Likelihood Ratio Test: \n")
    cat("Empirical Chi2 : \t", deparse(round(x$emp_Chi2,3)), "\n")
    cat("df: \t\t\t\t\t\t\t\t\t\t\t\t\t\t", x$df, "\n")
    cat("pvalue: \t\t\t\t\t\t\t\t\t\t", deparse(round(x$pval,3)), "\n")
  }  
  
}
