#' @title Set components names for vector names
#'
#' @description
#' If param has \code{cnames} set component names in a value.
#' Otherwise \code{x} is left unchanged.
#'
#' @template arg_par_or_set
#' @param x [any] \cr
#'   Param value(s).
#'   For a parameter set this must be a list in the correct order.
#' @return \code{x} with changed names.
#' @export
setValueCNames = function(par, x) {
  UseMethod("setValueCNames")
}

#' @export
setValueCNames.Param = function(par, x) {
  # do not set names for missing / req.params
  if (par$type %in% c("numericvector", "integervector", "logicalvector") && !isScalarNA(x))
    names(x) = par$cnames
  return(x)
}

#' @export
setValueCNames.ParamSet = function(par, x) {
  assertList(x, len = length(par$pars))
  Map(setValueCNames.Param, par$pars, x)
}
