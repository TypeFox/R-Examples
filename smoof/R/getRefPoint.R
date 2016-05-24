#' Returns the reference point of a multi-objective function.
#'
#' @template arg_smoof_function
#' @return [\code{numeric}]
#' @note Keep in mind, that this method makes sense only for multi-objective target functions.
#' @export
getRefPoint = function(fn) {
  UseMethod("getRefPoint")
}

#' @export
getRefPoint.smoof_single_objective_function = function(fn) {
  stopf("No reference point available for single-objective function.")
}

#' @export
getRefPoint.smoof_multi_objective_function = function(fn) {
  return(attr(fn, "ref.point"))
}

#' @export
getRefPoint.smoof_wrapped_function = function(fn) {
  return(getRefPoint(getWrappedFunction(fn)))
}
