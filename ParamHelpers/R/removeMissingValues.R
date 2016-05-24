#' Removes all scalar NAs from a parameter setting list.
#'
#' @param x [\code{list}]\cr
#'   List of paramter values.
#' @return [\code{list}].
#' @export
removeMissingValues = function(x) {
  return(Filter(Negate(isScalarNA), x))
}
