#' Repairs values of numeric and integer parameters out side of constraints.
#'
#' Clips values outside of box constraints to bounds.
#'
#' @template arg_parset
#' @param x [\code{list}]\cr
#'   List of parameter values. Must be in correct order.
#'   Values corresponding to non-numeric/integer types are left unchanged.
#' @param warn [\code{logical(1)}]\cr
#'   Boolean indicating whether a warning should be printed each time a value is repaired.
#'   Default is \code{FALSE}.
#' @return [\code{list}]:
#'   List of repaired points.
#' @export
repairPoint = function(par.set, x, warn = FALSE) {
  assertClass(par.set, "ParamSet")
  assertList(x)
  assertFlag(warn)
  Map(function(par, val) {
    if (isNumeric(par, include.int = TRUE)) {
      # repair scalar non-NAs and vectors
      if (!(length(val) == 1L && is.na(val)) && any(val < par$lower | val > par$upper)) {
        if (warn) {
          warningf("Repairing value for %s: %s", par$id, paramValueToString(par, val))
        }
        val = pmax(par$lower, val)
        val = pmin(par$upper, val)
      }
    }
    return(val)
  }, par.set$pars, x)
}
