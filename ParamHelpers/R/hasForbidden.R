#' Check parameter set for forbidden region.
#'
#' @template arg_parset
#' @return [\code{logical(1)}].
#' @export
hasForbidden = function(par.set) {
  assertClass(par.set, "ParamSet")
  return(!is.null(par.set$forbidden))
}


