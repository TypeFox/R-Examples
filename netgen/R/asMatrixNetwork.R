#' Convert network to matrix.
#'
#' @note If the instance contains of \eqn{n} depots, the depot coordinates fill the
#'   first \eqn{n} rows of the matrix.
#' @template arg_network
#' @param ... [any]\cr
#'   Currently not used.
#' @return [\code{matrix}]
#' @export
as.matrix.Network = function(x, ...) {
  if (!hasDepots(x)) {
    return(x$coordinates)
  }
  return(rbind(x$depot.coordinates, x$coordinates))
}
