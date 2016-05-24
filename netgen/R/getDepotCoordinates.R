#' Get coordinates of depots.
#'
#' @template arg_network
#' @return [\code{matrix}]
#' @export
getDepotCoordinates = function(x) {
  if (!hasDepots(x)) {
    stop("Object has no depots.")
  }
  x$depot.coordinates
}
