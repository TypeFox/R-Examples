#' Determines used factor levels.
#'
#' Determines the factor levels of a factor type vector
#' that are actually occuring in it.
#'
#' @param x [\code{factor}]\cr
#'   The factor.
#' @return [\code{character}]
#' @export
getUsedFactorLevels = function(x) {
  intersect(levels(x), unique(x))
}
