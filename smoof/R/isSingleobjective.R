#' Checks whether the given function is single-objective.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}] \code{TRUE} if function is single-objective.
#' @export
isSingleobjective = function(fn) {
  UseMethod("isSingleobjective")
}

#' @export
isSingleobjective.smoof_function = function(fn) {
  return(attr(fn, "n.objectives") == 1L)
}

#' @export
isSingleobjective.smoof_wrapped_function = function(fn) {
  return(isSingleobjective(getWrappedFunction(fn)))
}
