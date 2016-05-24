#' Helper function for determining the vector of attribute names
#' of a given object.
#'
#' @param obj [any]\cr
#'   Source object.
#' @return [\code{character}]
#'   Vector of attribute names for the source object.
#' @export
getAttributeNames = function(obj) {
  names(attributes(obj))
}
