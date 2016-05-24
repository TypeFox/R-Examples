#' Return the number of function evaluations performed by the wrapped
#' \code{smoof_function}.
#'
#' @param fn [\code{smoof_counting_function}]\cr
#'   Wrapped \code{smoof_function}.
#' @return [\code{integer(1)}]
#' @export
getNumberOfEvaluations = function(fn) {
  UseMethod("getNumberOfEvaluations")
}

#' @export
getNumberOfEvaluations.smoof_counting_function = function(fn) {
  return(environment(fn)$n.evals)
}

#' @export
getNumberOfEvaluations.smoof_wrapped_function = function(fn) {
  getNumberOfEvaluations(getWrappedFunction(fn))
}
