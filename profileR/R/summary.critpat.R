#' Summary of criterion-related profile analysis
#'
#' Provides a summary of the criterion-related profile analysis
#'
#' @method summary critpat
#' @param object object of class \code{critpat}
#' @param ... additional arguments affecting the summary produced.
#' @export
#' @seealso \code{\link{cpa}}
summary.critpat <- function(object, ...){
	cat("Call:\n")
	print(object$call)
  cat("\nRelability\n")
	print(object$r2)
  
  if(is.null(object$lvl.comp) == FALSE){
  cat("\n Level Component\n")
  print(object$lvl.comp)}
  if(is.null(object$pat.comp) == FALSE){
  cat("\n Pattern Component \n")
  print(object$pat.comp)
  }
}
