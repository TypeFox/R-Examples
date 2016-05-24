#' Checks whether the given function accept \dQuote{vectorized} input.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}]
#' @export
isVectorized = function(fn) {
  UseMethod("isVectorized")
}

#' @export
isVectorized.smoof_function = function(fn) {
  return(attr(fn, "vectorized"))
}

#' @export
isVectorized.smoof_wrapped_function = function(fn) {
  return(isVectorized(getWrappedFunction(fn)))
}
