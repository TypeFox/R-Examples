#' Get lower / upper bounds and allowed discrete values for parameters.
#'
#' \code{getLower} and \code{getUpper} return a numerical vector of lower and upper
#' bounds, \code{getValues} returns a list of possible value sets for discrete parameters.
#'
#' Parameters for which such bound make no sense - due to their type - are not present in the result.
#'
#' @template arg_parset
#' @param with.nr [\code{logical(1)}]\cr
#'   Should number from 1 to length be appended to names of vector params?
#'   Default is \code{FALSE}.
#' @return [\code{vector} | \code{list}]. Named by parameter ids.
#' @export
#' @examples
#' ps = makeParamSet(
#'   makeNumericParam("u"),
#'   makeIntegerParam("v", lower = 1, upper = 2),
#'   makeDiscreteParam("w", values = 1:2),
#'   makeNumericVectorParam("x", len = 2, lower = c(0, 10), upper = c(1, 11))
#' )
#' getLower(ps)
#' getUpper(ps)
#'
#' ps = makeParamSet(
#'   makeNumericParam("u"),
#'   makeDiscreteParam("v", values = c("a", "b")),
#'   makeDiscreteParam("w", values = list(a = list(), b = NULL))
#' )
#' getValues(ps)
getLower = function(par.set, with.nr = FALSE) {
  return(getBounds(par.set, type.of.bounds = "lower", with.nr = with.nr))
}

#' @export
#' @rdname getLower
getUpper = function(par.set, with.nr = FALSE) {
  return(getBounds(par.set, type.of.bounds = "upper", with.nr = with.nr))
}

#' @export
#' @rdname getLower
getValues = function(par.set) {
  assertClass(par.set, "ParamSet")
  types = getParamTypes(par.set)
  is.disc = types %in% c("discrete", "discretevector", "logical", "logicalvector")
  if (!any(is.disc))
    return(list())
  lapply(par.set$pars[is.disc], function(p) p$values)
}

# common functionality of getLower and getUpper
getBounds = function(par.set, type.of.bounds, with.nr = FALSE) {
  assertClass(par.set, "ParamSet")
  types = getParamTypes(par.set)
  is.num = types %in% c("numeric", "integer", "numericvector", "integervector")
  if (!any(is.num))
    return(numeric(0))
  bounds = lapply(par.set$pars[is.num], function(p) p[[type.of.bounds]])
  bounds = do.call(c, bounds)
  names(bounds) = getParamIds2(par.set$pars[is.num], repeated = TRUE, with.nr = with.nr)
  return(bounds)
}
