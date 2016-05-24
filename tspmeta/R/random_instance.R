#' Generates a random TSP instance by scattering random points in a hypercube.
#'
#' @param size [\code{integer(1)}]\cr
#'   Number of cities.
#' @param d [\code{integer(1)}]\cr
#'   Space dimensionality, e.g. 2D.
#'   Default is 2D.
#' @param lower [\code{numeric(1)}]\cr
#'   Lower box constraint for hypercube.
#'   Default is 0.
#' @param upper [\code{numeric(1)}]\cr
#'   upper box constraint for hypercube.
#'   Default is 1.
#' @return [\code{\link{tsp_instance}}].
#' @export
random_instance = function(size, d = 2, lower = 0, upper = 1) {
    x = runif(size * d, min = lower, max = upper)
	tsp_instance(coords = matrix(x, ncol = d))
}