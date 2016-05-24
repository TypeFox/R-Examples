#' Checks whether the objective function has constraints.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}]
#' @export
hasConstraints = function(fn) {
  UseMethod("hasConstraints")
}

#' @export
hasConstraints.smoof_function = function(fn) {
  (hasBoxConstraints(fn) || hasOtherConstraints(fn))
}

#' @export
hasConstraints.smoof_wrapped_function = function(fn) {
  fn = getWrappedFunction(fn)
  return(hasBoxConstraints(fn) || hasOtherConstraints(fn))
}

#' Checks whether the objective function has box constraints.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}]
#' @export
hasBoxConstraints = function(fn) {
  UseMethod("hasBoxConstraints")
}

#' @export
hasBoxConstraints.smoof_function = function(fn) {
  return(hasFiniteBoxConstraints(getParamSet(fn)))
}

#' @export
hasBoxConstraints.smoof_wrapped_function = function(fn) {
  return(hasFiniteBoxConstraints(getParamSet(getWrappedFunction(fn))))
}

#' Checks whether the objective function has other constraints.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}]
#' @export
hasOtherConstraints = function(fn) {
  UseMethod("hasOtherConstraints")
}

#' @export
hasOtherConstraints.smoof_function = function(fn) {
  return(!is.null(attr(fn, "constraint.fn")))
}

#' @export
hasOtherConstraints.smoof_wrapped_function = function(fn) {
  return(hasOtherConstraints(getWrappedFunction(fn)))
}
