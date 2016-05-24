#' Morphing of two networks with a convex combination of the coordinates.
#'
#' This function takes two (clustered) networks with equal number of nodes and,
#' if present, equal number of depots, and generates another instance by applying
#' a convex combination to the coordinates of node pairs. The node pairs are
#' determined by a point matching algorithm, which solves this assignement problem
#' via a integer programming procedure.
#' If both instances contain depots, point matching is done separately on depots
#' and the remaining nodes.
#'
#' @template arg_first_network
#' @template arg_second_network
#' @param alpha [\code{numeric(1)}]\cr
#'   Coeffiecient alpha for convex combination.
#' @param point.matching [\code{matrix} | NULL]\cr
#'   Point matching which shall be used for morphing. If \code{NULL}, an optimal
#'   point matching is generated via function \code{\link{getOptimalPointMatching}}.
#'   Default is \code{NULL}. Currently it is just possible to pass a point matching
#'   for instances without depots.
#' @param point.matching.algorithm [\code{function}]\cr
#'   Algorithm used to find a point matching. Default is \code{\link{getOptimalPointMatching}}.
#' @return [\code{Network}]
#'   Morphed network
#' @examples
#' x = generateRandomNetwork(n.points = 40L, n.depots = 2L)
#' y = generateClusteredNetwork(n.points = 40L, n.cluster = 2L, n.depots = 2L)
#' z = morphInstances(x, y, alpha = 0.2)
#' \dontrun{
#' library(gridExtra)
#' plot.list = list(autoplot(x), autoplot(z), autoplot(y))
#' plot.list$nrow = 1
#' do.call(grid.arrange, plot.list)
#' }
#' @seealso \code{\link{visualizeMorphing}}, \code{\link{visualizePointMatching}}
#' @export
morphInstances = function(x, y, alpha,
  point.matching = NULL,
  point.matching.algorithm = getOptimalPointMatching) {
  assertClass(x, "Network")
  assertClass(y, "Network")
  assertNumber(alpha, lower = 0, upper = 1, na.ok = FALSE)
  assertFunction(point.matching.algorithm)

  getPointMatchingAndMorphCoordinates = function(coords1, coords2) {
    point.matching = point.matching.algorithm(coords1, coords2)
    coordinates = makeConvexCombination(coords1, coords2[point.matching[, 2], ], alpha)
    return(coordinates)
  }

  x.coordinates = x$coordinates
  y.coordinates = y$coordinates

  if ((hasDepots(x) && !hasDepots(y)) || (!hasDepots(x) && hasDepots(y))) {
    stopf("Both or none of the instances must have depots")
  }

  if (hasDepots(x) && hasDepots(y) && !is.null(point.matching)) {
    stopf("Point matching parameter only supported for instances without depots!")
  }

  if (!is.null(point.matching)) {
    n = nrow(x.coordinates)
    assertMatrix(point.matching, ncols = 2L, nrows = n, any.missing = FALSE)
    if (any(point.matching[, 1] != seq(n))) {
      stopf("First column of 'point.matching' must contain the node IDs 1, ..., %i in this order.", n)
    }
    if (any(sort(point.matching[, 2]) != seq(n))) {
      stopf("Second column of 'point.matching' must contain a permutation of the node IDS 1, ..., %i.", n)
    }
    coordinates = makeConvexCombination(x.coordinates, y.coordinates[point.matching[, 2], ], alpha)
  } else {
    coordinates = getPointMatchingAndMorphCoordinates(x.coordinates, y.coordinates)
  }
  depot.coordinates = NULL

  if (all(hasDepots(x), hasDepots(y))) {
    x.n.depots = getNumberOfDepots(x)
    y.n.depots = getNumberOfDepots(y)
    if (x.n.depots != y.n.depots) {
      stopf("Number of depots must be equal, but x has %i and y has $i depots.",
        x.n.depots, y.n.depots)
    }
    x.depot.coordinates = getDepotCoordinates(x)
    y.depot.coordinates = getDepotCoordinates(y)
    depot.coordinates = getPointMatchingAndMorphCoordinates(x.depot.coordinates, y.depot.coordinates)
  }
  z = makeNetwork(
    coordinates = coordinates,
    depot.coordinates = depot.coordinates,
    lower = x$lower,
    upper = x$upper
  )

  attr(z, "morphed") = TRUE
  attr(z, "morphing.grade") = alpha
  return(z)
}
