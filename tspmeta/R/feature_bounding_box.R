#' Bounding box features.
#'
#' Determines the ratio of cities which lie within a certain
#' distance to the bounding box.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param distance_fraction [\code{numeric(1)}]\cr
#'   Distance ratio to bounding box.
#' @return [\code{list}].
#' @export
feature_bounding_box = function(x, distance_fraction = 0.1) {
	coords = x$coords
	ranges = apply(coords, 2, range)
	distances = diff(ranges)
	distance_x = distances[1] * distance_fraction
	distance_y = distances[2] * distance_fraction
	x_min = ranges[1,1] + distance_x
	x_max = ranges[2,1] - distance_x
	y_min = ranges[1,2] + distance_y
	y_max = ranges[2,2] - distance_y
	idx = which(coords[, 1] < x_min | coords[, 1] > x_max | coords[, 2] < y_min | coords[, 2] > y_max)
	if (length(idx) == 0) {
		idx = 0
	}
	res = list(ratio_of_cities_outside_box = length(unique(idx)) / number_of_cities(x))
	prefix = sprintf("bounding_box_%02i", floor(distance_fraction * 100))
	names(res) = paste(prefix, names(res), sep = "_")
	res
}
