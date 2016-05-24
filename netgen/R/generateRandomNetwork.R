#' Generates a random graph in a hypercube.
#'
#' @param n.points [\code{integer(1)}]\cr
#'   Number of points.
#' @param n.dim [\code{integer(1)}]\cr
#'   Number of dimensions. Default ist 2.
#' @param n.depots [\code{integer(1)}]\cr
#'   Number of depots in instances for the Vehicle Routing Problem (VRP).
#'   Default is NULL, i. e., no depots. The proceeding is as follows:
#'   If \code{n.depots} is 1, a random cluster center is defined to be the depot.
#'   If \code{n.depots} is 2, the second depot has maximal distance to the first.
#'   By convention the depots are placed as the first nodes in the coordinates
#'   matrix.
#' @param lower [\code{numeric(1)}]\cr
#'   Lower box constraint of cube.
#' @param upper [\code{numeric(1)}]\cr
#'   Upper box constraint of cube. Default is 100.
#' @template arg_name
#' @return [\code{Network}]
#' @examples
#' x = generateRandomNetwork(n.points = 100L, n.depots = 2L, upper = 50)
#' @export
generateRandomNetwork = function(n.points, n.dim = 2L, n.depots = NULL,
  lower = 0, upper = 100, name = NULL) {
  assertCount(n.points, na.ok = FALSE)
  assertInteger(n.dim, len = 1L, any.missing = FALSE, lower = 2L)

  if (!is.null(n.depots)) {
    assertInteger(n.depots, len = 1L, lower = 1L, upper = 2L)
  }

  assertNumber(lower, lower = 0, finite = TRUE)
  assertNumber(upper, finite = TRUE)
  if (!is.null(name)) assertCharacter(name, len = 1L, any.missing = FALSE)

  if (upper <= lower) {
    stopf("Argument 'upper' must be greater than argument 'lower'.")
  }

  coordinates = runif(n.points * n.dim, min = lower, max = upper)
  coordinates = matrix(coordinates, ncol = n.dim)

  depot.coordinates = NULL

  if (!is.null(n.depots)) {
    depot.coordinates = generateClusterCenters(n.cluster = 2L, lower = lower, upper = upper)
    if (n.depots == 1L) {
      depot.coordinates = depot.coordinates[1, , drop = FALSE]
    }
  }

  makeNetwork(
    name = coalesce(name, paste("RANDOM_", generateName(n.points, n.dim), sep = "")),
    coordinates = coordinates,
    depot.coordinates = depot.coordinates,
    lower = lower,
    upper = upper
  )
}
