#' Nearest neighbor features.
#'
#' Statistics describing the distribution of distances of each city
#' to its nearest neighbor.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{list}].
#' @export
feature_nnds = function(x) {
  d = x$dists
  diag(d) = Inf
  numvec_feature_statistics(apply(d, 1, min), "nnds")
}
