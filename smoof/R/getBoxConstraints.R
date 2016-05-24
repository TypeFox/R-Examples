#' Return lower box constaints.
#'
#' @template arg_smoof_function
#' @return [\code{numeric}]
#' @export
getLowerBoxConstraints = function(fn) {
  UseMethod("getLowerBoxConstraints")
}

#' @export
getLowerBoxConstraints.smoof_function = function(fn) {
  getLower(getParamSet(fn), with.nr = TRUE)
}

#' @export
getLowerBoxConstraints.smoof_wrapped_function = function(fn) {
  return(getLowerBoxConstraints(getWrappedFunction(fn)))
}

#' Return upper box constaints.
#'
#' @template arg_smoof_function
#' @return [\code{numeric}]
#' @export
getUpperBoxConstraints = function(fn) {
  UseMethod("getUpperBoxConstraints")
}

#' @export
getUpperBoxConstraints.smoof_function = function(fn) {
  getUpper(getParamSet(fn), with.nr = TRUE)
}

#' @export
getUpperBoxConstraints.smoof_wrapped_function = function(fn) {
  return(getUpperBoxConstraints(getWrappedFunction(fn)))
}
