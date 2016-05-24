#' Calculates list of all TSP features for an instance.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance
#' @param rescale [\code{logical(1)}]\cr
#'   Rescale \code{x} to \eqn{[0,1]^2} before calculation of features?
#'   Default is \code{TRUE}.
#' @return [\code{list}].
#'
#' @seealso \code{\link{feature_angle}},
#'          \code{\link{feature_centroid}},
#'          \code{\link{feature_cluster}},
#'          \code{\link{feature_bounding_box}},
#'          \code{\link{feature_chull}},
#'          \code{\link{feature_distance}},
#'          \code{\link{feature_modes}},
#'          \code{\link{feature_mst}},
#'          \code{\link{feature_nnds}}
#'
#' @export
#' @examples
#' x = random_instance(10)
#' print(features(x))
features = function(x, rescale = TRUE) {
	if (rescale) {
		x = rescale_instance(x)
	}
	c(feature_angle(x),
		feature_centroid(x),
		feature_cluster(x, 0.01),
		feature_cluster(x, 0.05),
		feature_cluster(x, 0.1),
		feature_bounding_box(x, 0.1),
		feature_bounding_box(x, 0.2),
		feature_bounding_box(x, 0.3),
		feature_chull(x),
		feature_distance(x),
		feature_modes(x),
		feature_mst(x),
		feature_nnds(x)
	)
}