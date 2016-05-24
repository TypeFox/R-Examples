##
## domination.R - Anything to do with Pareto dominance
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

##' \code{is_dominated} returns which points from a set are dominated
##' by another point in the set. \code{\%dominates\%} returns true if
##' \code{x} Pareto dominates \code{y} and
##' \code{is_maximally_dominated} returns TRUE for those points which
##' do not dominate any other points.
##'
##' @param points Matrix containing points one per column.
##'
##' @return For \code{is_dominated} and \code{is_maximally_dominated}
##' a boolean vector and for \code{\%dominates\%} a single boolean.
##'
##' @title Pareto dominance checks.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##'
##' @keywords optimize
##' @export
##' @rdname dom_op
is_dominated <- function(points) {
  #stopifnot(is.matrix(points))
  #stopifnot(is.numeric(points))
  .Call(do_is_dominated, points)
}

##' @export
##' @rdname dom_op
is_maximally_dominated <- function(points) {
  ## We should investiate a fast C implementation for this
  r <- nds_rank(points)
  r == max(r)
}

##' Nondominated points
##'
##' Return those points which are not dominated by another point in
##' \code{points}. This is the Pareto front approximation of the
##' point set. 
##'
##' @param points Matrix of points, one per column.
##' @return Those points in \code{points} which are not dominated by
##'   another point.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @export
##' @keywords optimize
nondominated_points <- function(points)
  points[,!is_dominated(points), drop=FALSE]

##' Nondominated sorting ranks
##'
##' Perform (partial) nondominated sort of the points in \code{points} and
##' return the rank of each point.
##'
##' @param points Matrix containing points one per column.
##' @param partial Optional integer specifying the number of points for
##'   which the rank should be calculated. Defaults to all points.
##'
##' @return Vector containing the ranks of the first \code{partial}
##'   individuals or all individuals.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @keywords optimize
##' @export
nds_rank <- function(points, partial) {
  #stopifnot(is.matrix(points))
  #stopifnot(is.numeric(points))
  
  if (missing(partial))
    partial <- ncol(points)
  else if (is.numeric(partial))
    partial <- as.integer(partial)
  else
    stopifnot(is.integer(partial))
  
  .Call(nondominated_order, points, partial)
}

##' @export
##' @rdname nds_rank
nondominated_ordering <- function(points, partial) {
  .Deprecated("nds_rank")
  nds_rank(par, partial)
}
