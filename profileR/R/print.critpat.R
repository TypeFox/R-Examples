#' Print a criterion-related profile analysis
#' 
#' Prints the default output from fitting the \code{cpa} function.
#' @param x object of class \code{critpat} returned from the \code{cpa} function
#' @param ... additional objects of the same type.
#' @method print critpat
#' @seealso \code{\link{cpa}}
#' @export

print.critpat <- function(x, ...){
  if (!inherits(x, "critpat")) stop("Use only with 'critpat' objects.\n")
	cat("Call:\n")
	print(x$call)
  cat("\nCoefficients\n")
	print(x$b)
  invisible(x) 
}
