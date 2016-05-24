#' Check if network has depots.
#'
#' @template arg_network
#' @return [\code{logical(1)}]
#' @export
hasDepots = function(x) {
  assertClass(x, "Network")
  !is.null(x$depot.coordinates)
}
