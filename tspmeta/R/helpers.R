# L2-Norm.
#
# @param x[\code{numeric}]\cr
#   Numeric vector.
#
# @return [\code{numeric(1)}].
l2_norm = function(x) {
	sqrt(sum(x * x))
}

# Euclidean distance.
#
# @param x [\code{numeric}]\cr
#   Coordinates of a city.
# @param y [\code{matrix}]\cr
#   Numeric matrix of city coordinates,
#   rows denote cities.
#
# @return [\code{numeric}]
#   Euclidean distances between city and all other cities.
eucl_distance = function(x, y){
	sapply(1:nrow(y), function(i) {
		sqrt(crossprod(x - y[i, ]))
	})
}

# Angle between the segments 'ab' and 'ac'.
#
# @param a [\code{numeric}]\cr
#   Central point.
# @param b [\code{numeric}]\cr
#   First neighbor.
# @param c [\code{numeric}]\cr
#   Second neighbor.
#
# @return [\code{numeric(1)}]
#   The (smaller) unsigned angle between 'ab' and 'ac'.
angle_between_points = function(a, b, c) {
	v1 = (a - b) / l2_norm(a - b)
	v2 = (a - c) / l2_norm(a - c)
	z = sum(v1 * v2)
	# FIXME: check this again!
	if (1 - z < 1e-10)
		angle = 0
	else if	(z - (-1) < 1e-10)
		angle = pi
	else
		angle = acos(z)
	#angle = acos(sum(v1 * v2))
	min(abs(angle + c(0, 2 * pi, -2 * pi)))
}

#' Return the center of all cities of a TSP instance.
#'
#' @param instance [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#'
#' @return [\code{numeric(2)}] Center of all cities of the TSP instance.
#'
#' @export
center_of_mass = function(instance) {
	assertClass(instance, "tsp_instance")
	colMeans(instance$coords)
}
