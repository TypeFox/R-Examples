#' Angle features.
#'
#' Statistics of the distribution of the angle between a node and its 2 next neighbors.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{list}].
#' @export
#FIXME: only works for 2D
feature_angle = function(x) {
  #FIXME: delete duplicate cities; otherwise computation of angles fails
  coords = x$coords
  angles = sapply(1:number_of_cities(x), function(city) {
	neighbor = order(x$dists[city,])[2:3]  # 2 next neighbors.
	angle_between_points(coords[city, ], coords[neighbor[1], ],
	  coords[neighbor[2], ])
  })
	numvec_feature_statistics(angles, "angle")
}
