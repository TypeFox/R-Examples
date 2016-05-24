#' Summarize an epplab Object
#' 
#' Summarizes and prints an \code{epplab} object in an informative way.
#' 
#' The option \code{which} can restrict the output to certain simulation runs.
#' In case of many simulations, this might improve the readability.
#' 
#' @name summary.epplab
#' @aliases summary.epplab summary-method summary,epplab-method
#' @docType methods
#' @param object Object of class \code{epplab}.
#' @param which Summary for \code{which} simulation runs
#' @param ... Additional parameters
#' @author Daniel Fischer
#' @keywords methods print
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMin",n.simu=10, maxiter=20)
#' summary(res)
#' 
#' @export
summary.epplab <- function(object, which=1:10, ...){
  which <- which[which<=length(object$PPindexVal)]
  cat("REPPlab Summary\n")
  cat("---------------\n")
  cat("Index name       :",object$PPindex,"\n")
  cat("Index values     :",object$PPindexVal[which],"\n")
  cat("Algorithm used   :",object$PPalg,"\n")
  cat("Sphered          :",object$sphered,"\n")
  cat("Iterations       :",object$PPiter[which],"\n")
  invisible(object)
} 
