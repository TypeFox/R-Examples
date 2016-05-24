##' Determine which points are on the edge of a Pareto-front approximation.
##'
##' @param front Pareto-front approximation.
##' 
##' @return An integer vector containing the indicies of the points
##' (columns) of \code{front} which are on the edge of the
##' Pareto-front approximation.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##'
##' @export
which_points_on_edge <- function(front) {
  which(.Call(do_which_points_on_edge, front))
}
