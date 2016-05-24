#' Checks whether the given function is noisy.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}]
#' @export
isNoisy = function(fn) {
  UseMethod("isNoisy")
}

#' @export
isNoisy.smoof_function = function(fn) {
  return(attr(fn, "noisy"))
}

#' @export
isNoisy.smoof_wrapped_function = function(fn) {
  return(isNoisy(getWrappedFunction(fn)))
}
