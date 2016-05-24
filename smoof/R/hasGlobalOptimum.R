#' Checks whether global optimum is known.
#'
#' @template arg_smoof_function
#' @return [\code{logical(1)}]
#' @export
hasGlobalOptimum = function(fn) {
  UseMethod("hasGlobalOptimum")
}

#' @export
hasGlobalOptimum.smoof_single_objective_function = function(fn) {
  return(!is.null(attr(fn, "global.opt.params")))
}

#' @export
hasGlobalOptimum.smoof_multi_objective_function = function(fn) {
  return(FALSE)
}

#' @export
hasGlobalOptimum.smoof_wrapped_function = function(fn) {
  return(hasGlobalOptimum(getWrappedFunction(fn)))
}
