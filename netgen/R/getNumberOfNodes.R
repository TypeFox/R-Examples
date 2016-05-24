#' Returns number of nodes of a network.
#'
#' @template arg_network
#' @return [\code{integer(1)}]
#'   Number of nodes of the network.
#' @export
getNumberOfNodes = function(x) {
  assertClass(x, "Network")
  nrow(x$coordinates)
}
