#' @title Computes optimal point assignment for two sets of points of equal size.
#'
#' @description Internally it handles the points and the possible matchings as a bi-partite
#' graphs and finds an optimal matching due to euclidean distance by an
#' efficient linear programming solver.
#'
#' @param x [\code{Network} | \code{matrix}]\cr
#'   First network or matrix of coordinates of the first point set.
#' @param y [\code{Network} | \code{matrix}]\cr
#'   Second network or matrix of coordinates of the second point set.
#' @param method [\code{character(1)}]\cr
#'   Method used to solve the assignment problem. There are currently two methods
#'   available:
#'   \describe{
#'     \item{lp}{Solves the problem be means of linear programming with the
#'     \pkg{lpSolve} package. This is the default.}
#'     \item{push_relabel}{The assignment problem can be formulated as a
#'     matching problem on bipartite graphs. This method makes use of the
#'     push-relabel algorithm from the \pkg{igraph}.}
#'     \item{random}{Random point matching.}
#'   }
#' @return [\code{matrix}]
#'   Each row consists of the indizes of the pairwise matchings.
#' @seealso \code{\link{visualizePointMatching}}
#' @export
getOptimalPointMatching = function(x, y, method = "lp") {
  coords1 = x
  coords2 = y
  if (isNetwork(x)) {
    coords1 = x$coordinates
  }
  if (isNetwork(y)) {
    coords2 = y$coordinates
  }
  assertMatrix(coords1, mode = "numeric")
  assertMatrix(coords2, mode = "numeric")
  if (ncol(coords1) > 2L || ncol(coords2) > 2L) {
    stopf("Point matching: At the moment only 2-dimensional point sets can be matched.")
  }

  if (any(dim(coords1) != dim(coords2))) {
    stopf("Point matching: Both coordinate matrizes need to have the same dimension.")
  }

  mapping = list(
    "lp" = getPointMatchingBySolvingLP,
    "push_relabel" = getPointMatchingByPushRelabelAlgorithm,
    "random" = getRandomPointMatching
  )

  assertChoice(method, choices = names(mapping))

  matching.algorithm = mapping[[method]]
  return(matching.algorithm(coords1, coords2))
}

# Solve assignement problem by means of linear programming with the lpSolve
# package.
getPointMatchingBySolvingLP = function(coords1, coords2) {
  dist.matrix = matrix(nrow = nrow(coords1), ncol = nrow(coords2))
  for (i in seq(nrow(coords1))) {
    for (j in seq(nrow(coords2))) {
      dist.matrix[i, j] = euklideanDistance(coords1[i, ], coords2[j, ])
    }
  }

  requirePackages("lpSolve", why = "netgen::getPointMatchingBySolvingLP")
  lp.res = lp.assign(dist.matrix)
  if (lp.res$status != 0) {
    stop("Failed to find LP solution! No point matching possible.")
  }
  lp.res = lp.res$solution

  # now construct mapping matrix
  res = matrix(nrow = nrow(lp.res), ncol = 2)
  res[, 1] = 1:nrow(lp.res)
  for (i in 1:nrow(lp.res)) {
    res[i, 2] = which(lp.res[i, ] != 0)
  }
  #FIXME: what the fuck is going on here??? Each line consists of exactly one 1
  # but R does not find it! if I ask for != 0, it works! -.-w
  #res[, 2] = as.numeric(apply(lp.res, 1, function(row) as.numeric(which(row == 1))))
  return(res)
}

# Uses push-relabel algorithm to comute a minimum weight matching on bipartite
# graph with the igraph package.
getPointMatchingByPushRelabelAlgorithm = function(coords1, coords2) {
  requirePackages("igraph", why = "netgen::getPointMatchingByPushRelabelAlgorithm")
  n = nrow(coords1)

  # generate a grid of paired node IDs (of bipartite graph)
  grid = expand.grid("x" = 1:n, "y" = 1:n)

  # add weight property, i.e., edge costs
  grid$weight = apply(grid, 1, function(row) {
    x.id = row[1]
    y.id = row[2]
    x.coord = coords1[x.id, ]
    y.coord = coords2[y.id, ]
    sqrt(sum((x.coord - y.coord)^2))
  })

  # since we want a minimum cost matching we need to adapt the weight
  grid$weight = max(grid$weight) - grid$weight

  # moreover we need distinct nodes sets (we add n here and substract it later)
  grid$y = grid$y + n

  # generate igraph bipartite graph
  gr = igraph::graph.data.frame(d = grid)

  # make graph bipartite
  V(gr)$type = rep(c(TRUE, FALSE), each = n)

  # compute matching
  matching = igraph::maximum.bipartite.matching(gr)

  # only the first part is interesting for us.
  matching = matching$matching[1:n]
  matching = matrix(c(as.integer(names(matching)), as.integer(matching)), ncol = 2)

  # revert "grid$y = grid$y + n"
  matching[, 2] = matching[, 2] - n
  return(matching)
}

# Simple random point assignment.
getRandomPointMatching = function(coords1, coords2) {
  ids = 1:nrow(coords1)
  matching = matrix(c(ids, sample(ids)), ncol = 2L)
  return(matching)
}
