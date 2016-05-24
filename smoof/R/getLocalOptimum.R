#' Returns the local optima.
#'
#' @description This function returns the parameters and objective values of
#' all local optima (including the global one).
#' @template arg_smoof_function
#' @return [\code{list}] List containing the following entries:
#' \itemize{
#'   \item{param [\code{list}]}{List of parameter values per local optima.}
#'   \item{value [\code{list}]}{List of objective values per local optima.}
#'   \item{is.minimum [\code{logical(1)}]}{Are the local optima minima or maxima?}
#' }
#' @note Keep in mind, that this method makes sense only for single-objective target functions.
#' @export
getLocalOptimum = function(fn) {
  UseMethod("getLocalOptimum")
}

#' @export
getLocalOptimum.smoof_single_objective_function = function(fn) {
  return(list(
    params = attr(fn, "local.opt.params"),
    values = attr(fn, "local.opt.value"),
    is.minimum = attr(fn, "minimize")
  ))
}

#' @export
getLocalOptimum.smoof_multi_objective_function = function(fn) {
  stopf("No local optimum available for multi-objective function.")
}

#' @export
getLocalOptimum.smoof_wrapped_function = function(fn) {
  return(getLocalOptimum(getWrappedFunction(fn)))
}
