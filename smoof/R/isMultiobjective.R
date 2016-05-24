#' Checks whether the given function is multi-objective.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}] \code{TRUE} if function is multi-objective.
#' @export
isMultiobjective = function(fn) {
  UseMethod("isMultiobjective")
}

#' @export
isMultiobjective.smoof_function = function(fn) {
  return(attr(fn, "n.objectives") >= 2L)
}

#' @export
isMultiobjective.smoof_wrapped_function = function(fn) {
  return(isMultiobjective(getWrappedFunction(fn)))
}
