#' Get parameter set.
#'
#' Each smoof function contains a parameter set of type \code{\link[ParamHelpers]{ParamSet}}
#' assigned to it, which describes types and bounds of the function parameters.
#' This function returns the parameter set.
#'
#' @template arg_smoof_function
#' @return [\code{\link[ParamHelpers]{ParamSet}}]
#' @examples
#' fn = makeSphereFunction(3L)
#' ps = smoof::getParamSet(fn)
#' print(ps)
#' @export
getParamSet = function(fn) {
  UseMethod("getParamSet")
}

#' @export
getParamSet.smoof_function = function(fn) {
  return(attr(fn, "par.set"))
}

#' @export
getParamSet.smoof_wrapped_function = function(fn) {
  return(getParamSet(getWrappedFunction(fn)))
}
