#' @rdname mprm
#'
#' @export
#'
#' @method summary MPRM
#'
#' @param object object of class \code{MPRM}
#' @param \dots \dots{}


summary.MPRM <-
function(object, ...){
  
  cat("\n Call: ", deparse(object$call), "\n\n")
  
  cat("Function calls: \n")
    print(object$fun_calls)
  
  cat("Convergence: \t", ifelse(object$convergence==0, "convergence", "no convergence"), "\n\n")
  
  cat("Deviance: \t", deparse(round(object$logLikelihood*(-2),3)), "\n")
  cat("Number of Parameters: \t", length(object$estpar), "\n\n")
  
  cat("-----------------------------------------------------\n")
  cat("Parameter estimates: \n")
  print(object$itempar)
  print(object$itempar_se)
  
}
