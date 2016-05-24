#' Generates a grid network.
#'
#' @param n.points.per.dim [\code{integer(1)}]\cr
#'   Number of points in each dimension.
#' @template arg_n_dim
#' @template arg_lower
#' @template arg_upper
#' @template arg_name
#' @return [\code{Network}]
#' @examples
#' x = generateGridNetwork(n.points.per.dim = 10L, upper = 50)
#' @note Grid networks with depots are not supported at the moment.
#' @export
generateGridNetwork = function(n.points.per.dim = NULL, n.dim = 2L, lower = 0, upper = 100, name = NULL) {
  assertCount(n.points.per.dim, na.ok = FALSE)
  assertInteger(n.dim, len = 1L, any.missing = FALSE, lower = 2L)
  assertNumber(lower, lower = 0, finite = TRUE)
  assertNumber(upper, finite = TRUE)

  if (!is.null(name)) {
    assertCharacter(name, len = 1L, any.missing = FALSE)
  }

  if (upper <= lower) {
    stopf("Argument 'upper' must be greater than argument 'lower'.")
  }

    # build source sequence for each dimension
  source.seq = seq(lower, upper, length.out = round(n.points.per.dim, digits = 0L))
  x = list()
  for (i in 1:n.dim) {
    x[[i]] = source.seq
  }

  # generate coordinates
  coordinates = as.matrix(do.call(expand.grid, x))

  makeNetwork(
    name = coalesce(name, paste("GRID_", generateName(n.points.per.dim ** 2, n.dim), sep = "")),
    coordinates = coordinates,
    lower = lower,
    upper = upper
    )
}
