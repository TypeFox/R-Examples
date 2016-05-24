#' Reset evaluation counter.
#'
#' @param fn [\code{smoof_counting_function}]\cr
#'   Wrapped \code{smoof_function}.
#' @export
resetEvaluationCounter = function(fn) {
  UseMethod("resetEvaluationCounter")
}

#' @export
resetEvaluationCounter.smoof_counting_function = function(fn) {
  environment(fn)$n.evals = 0L
}

#' @export
resetEvaluationCounter.smoof_wrapped_function = function(fn) {
  return(resetEvaluationCounter(getWrappedFunction(fn)))
}
