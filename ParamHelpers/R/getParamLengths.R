#' Return lengths of parameters in parameter set.
#'
#' Useful for vector parameters.
#'
#' @template arg_parset
#' @return [\code{integer}]. Named and in same order as \code{par.set}.
#' @examples
#' ps = makeParamSet(
#'   makeNumericParam("u"),
#'   makeIntegerParam("v", lower = 1, upper = 2),
#'   makeDiscreteParam("w", values = 1:2),
#'   makeDiscreteVectorParam("x", len = 2, values = c("a", "b"))
#' )
#' getParamLengths(ps)
#' # the length of the vector x is 2, for all other single value parameters the length is 1.
#' @export
getParamLengths = function(par.set) {
  assertClass(par.set, "ParamSet")
  # if we dont do this check we get an empty list
  if (isEmpty(par.set))
    return(integer(0))
  else
    return(extractSubList(par.set$pars, "len"))
}
