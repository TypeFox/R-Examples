#' A sub-class of the class \code{\linkS4class{mtkEvaluatorResult}} used to hold the results of the model simulation.
#'
#' @title The mtkSystemEvaluatorResult class
#' @exportClass mtkSystemEvaluatorResult

setClass("mtkSystemEvaluatorResult",
		
		contains=c("mtkEvaluatorResult")
)

#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Evaluator.
#'  @param information a named list to provide supplementary information about the model simulation  and its results.

#' @return an object of class \code{\linkS4class{mtkSystemEvaluatorResult}}
#' @examples mtkmtkSystemEvaluatorResult()
#' @export mtkmtkSystemEvaluatorResult
#' @title The constructor

mtkSystemEvaluatorResult <- function(main, information=NULL) {
	res <- new("mtkSystemEvaluatorResult", main=main, information=information)
	return(res)
}
