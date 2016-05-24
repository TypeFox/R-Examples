#' Anova Tables
#'
#' Computes an analysis of variance table for a criterion-related profile analysis
#'
#' @importFrom stats anova
#' @param object an object containing the results returned by a model fitting \code{cpa}.
#' @param ... additional objects of the same type.
#' @method anova critpat
#' @export
#' @seealso \code{\link{cpa}}
anova.critpat <- function(object, ...){
	cat("Call:\n")
	print(object$call)
	cat("\nAnalysis of Variance Table\n")
	cat("\n")
	print(object$ftable)
}
