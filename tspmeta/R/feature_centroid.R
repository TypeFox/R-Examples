#' Centroid features.
#' 
#' Includes the coordinates of the mean coordinates of the the point cloud
#' and the statistics of the distances of all cities from it.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{list}].
#' @export
feature_centroid = function(x) {
  coords = x$coords
  centroid = colMeans(coords)
  c_dists = apply(coords, 1, function(point) l2_norm(point - centroid))
  res = list(centroid_centroid_x = centroid[1], centroid_centroid_y = centroid[2])
  c(res, numvec_feature_statistics(c_dists, "centroid_dist"))
}
