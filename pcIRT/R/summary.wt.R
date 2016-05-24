#' @rdname wt
#'
#' @export
#'
#' @method summary wt
#'
#' @param object object of class \code{wt}
#' @param \dots \dots{}

summary.wt <-
function(object, ...){
  
  cat("\n Weight Test unidimensional polytomous model, \n\n")
  
  cat("Unconstrained model: \n")
  cat("LogLikelihood: ", object$unconstrLoglikelihood, "\n")
  cat("Number of parameters: ", object$unconstrNrPar, "\n\n")
  
  cat("Score parameters: \n")
  print(object$unconstrScoreParameter)
  
    cat("-----------------------------------------------------\n")
  
    cat("Constrained model: \n")
    cat("LogLikelihood: ", object$constrLogLikelihood, "\n")
    cat("Number of parameters: ", object$constrNrPar, "\n\n")
  
    cat("-----------------------------------------------------\n")
  
    cat("Likelihood Ratio Test: \n")
    cat("Empirical Chi2 : \t", deparse(round(object$emp_Chi2,3)), "\n")
    cat("df: \t\t\t\t\t\t\t\t\t\t\t\t\t\t", object$df, "\n")
    cat("pvalue: \t\t\t\t\t\t\t\t\t\t", deparse(round(object$pval,3)), "\n")
  
}
