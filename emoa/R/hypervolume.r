##
## hypervolume.r - Functions for calculating the dominated hypervolume
##
## Authors:
##   Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

##' Dominated Hypervolume calculation
##' 
##' \code{dominated_hypervolume} calculates the dominated hypervolume of
##' the points in \code{points}. 
##' 
##' \code{hypervolume_contribution} calculates the hypervolume
##' contribution of each point.
##' 
##' If no reference point \code{ref} is given, one is automatically
##' calculated by determening the maximum in each coordinate.
##'  
##' Currently only one general algorithm is implemented due to Fonseca
##' et.al. but work is underway to include others such as the Beume &
##' Rudolph approach as well as the approach by Bradstreet et.al.
##'  
##' The 1D and 2D cases are handle seperately by efficient algorithms.
##' Calculates the exact dominated hypervolume of the points given in
##' \code{x} subject to the reference point \code{ref}.
##' 
##' @param points Matrix containing the points one per column.
##' @param ref Optional reference point. If not provided the maximum
##'   in each dimension is used.
##' @return For \code{dominated_hypervolume} the dominated hypervolume
##'   by the points in \code{points} with respect to the reference point
##'   \code{ref}. For \code{hypervolume_contribution} a vector giving
##'   the hypervolume soley dominated by that point.
##' 
##' @seealso \code{\link{nondominated_points}} to extract the pareto
##' front approximation from a given set of points and
##' \code{\link{nds_hv_selection}} for a selection strategy based on
##' the hypervolume contribution of each point.
##'
##' @references
##' This code uses version 1.3 of the hypervolume code available from
##' \url{http://iridia.ulb.ac.be/~manuel/hypervolume}. For a
##' description of the algorithm see
##'
##' Carlos M. Fonseca, Luis Paquete, and Manuel Lopez-Ibanez. \emph{An
##' improved dimension-sweep algorithm for the hypervolume
##' indicator}. In IEEE Congress on Evolutionary Computation, pages
##' 1157-1163, Vancouver, Canada, July 2006.
##' 
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' 
##' @export
##' @keywords optimize
dominated_hypervolume <- function(points, ref) {
  ## Possibly infer reference point:
  if (missing(ref))
    ref <- apply(points, 1, max)

  ## Sanity checks:
  if (!is.matrix(points))
    stop("Pareto front must be a matrix")
  if (nrow(points) != length(ref))
    stop("Reference point and front must have the same dimension.")
  if (!all(is.finite(ref))) {
    warning("Reference point containes non finite values.")
    return(NaN)
  }
  if (!all(is.finite(points))) {
    warning("Front includes non finite points.")
    return(NaN)
  }
  
  .Call(do_dominated_hypervolume, points, ref, PACKAGE="emoa")
}

##' @export
##' @rdname dominated_hypervolume
hypervolume_contribution <- function(points, ref) {
  ## Possibly infer reference point:
  if (missing(ref))
    ref <- apply(points, 1, max) + 1

  ## Sanity checks:
  if (!is.matrix(points))
    stop("Pareto front must be a matrix")
  if (nrow(points) != length(ref))
    stop("Reference point and front must have the same dimension.")

  ## Call C code:
  .Call(do_hv_contrib, points, ref, PACKAGE="emoa")
}
