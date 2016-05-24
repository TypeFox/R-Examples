#' Get the number of clusters of a network.
#'
#' @template arg_network
#' @return [\code{integer(1)}]
#'   Number of clusters.
#' @note For simple random or grid networks this function always returns 1.
#' @export
getNumberOfClusters = function(x) {
  assertClass(x, "Network")
  n.cluster = 1L
  if (testClass(x, "ClusteredNetwork")) {
    # NOTE: depots are encoded as special cluster members of the depot 0.
    # we filter these out!
    membership = x$membership[which(x$membership != 0)]
    n.cluster = length(unique(membership))
  }
  n.cluster
}
