#' Returns the global optimum and its value.
#'
#' @template arg_smoof_function
#' @return [\code{list}] List containing the following entries:
#' \itemize{
#'   \item{param [\code{list}]}{Named list of parameter value(s).}
#'   \item{value [\code{numeric(1)}]}{Optimal value.}
#'   \item{is.minimum [\code{logical(1)}]}{Is the global optimum a minimum or maximum?}
#' }
#' @note Keep in mind, that this method makes sense only for single-objective target function.
#' @export
getGlobalOptimum = function(fn) {
  UseMethod("getGlobalOptimum")
}

#' @export
getGlobalOptimum.smoof_single_objective_function = function(fn) {
  return(list(
    param = attr(fn, "global.opt.params"),
    value = attr(fn, "global.opt.value"),
    is.minimum = attr(fn, "minimize")
  ))
}

#' @export
getGlobalOptimum.smoof_multi_objective_function = function(fn) {
  stopf("No global optimum available for multi-objective function.")
}

#' @export
getGlobalOptimum.smoof_wrapped_function = function(fn) {
  return(getGlobalOptimum(getWrappedFunction(fn)))
}
