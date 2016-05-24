#' Check whether the function is counting its function evaluations.
#'
#' In this case the function is of type \code{smoof_counting_function} or it
#' is further wrapped by another wrapper. This function then checks recursively,
#' if there is a counting wrapper.
#'
#' @param object [any]\cr
#'   Arbitrary R object.
#' @return \code{logical(1)}
#' @seealso \code{\link{addCountingWrapper}}
#' @export
doesCountEvaluations = function(object) {
  res = inherits(object, what = "smoof_counting_function")
  if (!res && isWrappedSmoofFunction(object)) {
    res = (res || doesCountEvaluations(getWrappedFunction(object)))
  }
  return(res)
}
