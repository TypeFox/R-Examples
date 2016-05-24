#' Determine the number of parameters.
#'
#' @template arg_smoof_function
#' @return [\code{integer(1)}]
#' @export
getNumberOfParameters = function(fn) {
  UseMethod("getNumberOfParameters")
}

#' @export
getNumberOfParameters.smoof_function = function(fn) {
  return(sum(getParamLengths(getParamSet(fn))))
}

#' @export
getNumberOfParameters.smoof_wrapped_function = function(fn) {
  return(getNumberOfParameters(getWrappedFunction(fn)))
}
