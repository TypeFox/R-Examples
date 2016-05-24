#' Returns information on the number of parameters of a each type.
#'
#' @template arg_parset
#' @return [\code{list}]
#'  Named list which contains for each supported parameter type the
#'  number of parameters of this type in the given ParamSet.
#' @export
getParamTypeCounts = function(par.set) {
  assertClass(par.set, "ParamSet")
  supported.types = getSupportedParamTypes()
  par.types = getParamTypes(par.set)
  count = lapply(supported.types, function(type) {
    sum(par.types == type)
  })
  names(count) = supported.types
  return(count)
}

# Returns a vector of supported parameter types.
getSupportedParamTypes = function() {
  return(
    c("numeric", "numericvector",
      "integer", "integervector",
      "discrete", "discretevector",
      "logical", "logicalvector",
      "character", "charactervector",
      "function",
      "untyped"))
}
