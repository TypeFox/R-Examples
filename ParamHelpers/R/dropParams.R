#' Drop Params from ParamSet by ids.
#'
#' @template arg_parset
#' @param drop [\code{character}]\cr
#'   \code{id}s of the \code{\link{Param}}s in the \code{\link{ParamSet}} to drop from the ParamSet.
#' @return [\code{\link{ParamSet}}].
#' @export
dropParams = function(par.set, drop) {
  assertClass(par.set, "ParamSet")
  assertSubset(drop, getParamIds(par.set))
  par.set$pars = Filter(function(p) p$id %nin% drop, par.set$pars)
  return(par.set)
}
