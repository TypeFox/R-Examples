#' Return number of parameters in set.
#'
#' Either number of parameters or sum over parameter lengths.
#'
#' @template arg_parset
#' @param devectorize [\code{logical(1)}]\cr
#'   Sum over length of vector parameters?
#'   Default is code{FALSE}.
#' @return [\code{integer}].
#' @examples
#' ps = makeParamSet(
#'   makeNumericParam("u"),
#'   makeDiscreteVectorParam("x", len = 2, values = c("a", "b"))
#' )
#' getParamNr(ps)
#' getParamNr(ps, devectorize = TRUE)
#' @export
getParamNr = function(par.set, devectorize = FALSE) {
  assertClass(par.set, "ParamSet")
  assertFlag(devectorize)
  if (devectorize)
    return(sum(getParamLengths(par.set)))
  else
    return(length(par.set$pars))
}
