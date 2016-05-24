#' Returns the number of depots of a network.
#'
#' @template arg_network
#' @return [\code{integer(1)}]
#' @export
getNumberOfDepots = function(x) {
  assertClass(x, "Network")
  if (hasDepots(x)) {
    return(nrow(x$depot.coordinates))
  }
  return(0L)
}
