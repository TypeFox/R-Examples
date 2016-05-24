#' @title Check if network is euclidean.
#'
#' @description Check if a \code{Network} object has euclidean coordinates.
#'
#' @template arg_network
#' @return [\code{logical(1)}]
#' @export
isEuclidean = function(x) {
  return(x$edge.weight.type %in% c("EUC_2D", "EUC_3D"))
}
