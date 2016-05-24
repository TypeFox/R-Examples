#' Greedy point matching
#'
#' Pairs of cities are matched in a greedy fashion for morphing,
#' first the closest pair w.r.t. euclidean distance, then the clostest
#' pair of the remaining cities, and so on.
#'
#' @param x [\code{tsp_instance}]\cr
#'   First TSP instance.
#' @param y [\code{tsp_instance}]\cr
#'   Second TSP instance.
#' @return [\code{matrix}]
#'   Numeric matrix of point indices with shortest distance.
#' @export
greedy_point_matching = function(x, y) {
  assertClass(x, "tsp_instance")
  assertClass(y, "tsp_instance")
  x = x$coords
  y = y$coords
  # all distances to y_i's from each point in x
  lists_dists = lapply(1:nrow(x), function(i) eucl_distance(x[i, ], y))
  len = length(lists_dists)

  # initialize index matrix
  ind = matrix(c(1:len, rep(0, len)), byrow = FALSE, ncol = 2, nrow = len)
  for (i in c(1:len)){
    ## index of point in x for which the minimum distance to all
    ## points in y is smallest
    m  = which.min(do.call(rbind,
                            lapply(lists_dists,
                                   function(x) min(x, na.rm = TRUE))))
                                        # index of point in y with smallest distance
    ind[m, 2] = which.min(lists_dists[[m]])

    ## discard the selected point in y
    for (j in c(1:len)[-m]){
      if((is.na(lists_dists[[j]][1]) == TRUE) | (is.finite(lists_dists[[j]][1]))){
        lists_dists[[j]][ind[m,2]] = NA
      }
    }

    ## set entry to Infinity to avoid repeated selection
    lists_dists[m] = Inf
  }
  ind
}

# Random point matching.
#
# @param x [\code{tsp_instance}]
#   First TSP instance.
# @param y [\code{tsp_instance}]
#   Second TSP instance.
# @param ntries [\code{integer(1)}] \cr
#   Number of random matchings to try. The best one is returned.
#   Default is 100.
# @return [\code{matrix}]
#   A matrix with 2 columns that specifies the best of the
#   \code{ntries} matching that was found.
random_point_matching = function(x, y, ntries = 100) {
  assertClass(x, "tsp_instance")
  assertClass(y, "tsp_instance")
  assertCount(x, positive = TRUE, na.ok = FALSE)

  best_matching = NULL
  best_matching_error = Inf
  n = nrow(x$coords)
  for (i in 1:ntries) {
    matching = cbind(1:n, sample(1:n))
    error = matching_error(matching, x, y)
    if (error < best_matching_error) {
      best_matching = matching
      best_matching_error = error
    }
  }
  best_matching
}

# Quantify error of instance matching.
#
# Uses the sum of the squared distance between the matched cities as
# a measure for the error of the matching. If a matching is perfect,
# the error would be zero, else greater than zero.
#
# @param matching [\code{matrix}]\cr
#   A matrix with 2 columns that matches cities of \code{x} to cities in \code{y}.
# @param x [\code{tsp_instance}]
#   First TSP instance.
# @param y [\code{tsp_instance}]
#   Second TSP instance.
# @return The error of the matching.
matching_error = function(matching, x, y) {
  D = x$coords[matching[, 1], ] - y$coords[matching[, 2], ]
  sum(D * D)
}
