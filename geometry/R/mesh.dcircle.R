##' Signed distance from points \code{p} to boundary of circle to
##' allow easy definition of regions in \code{\link{distmesh2d}}.
##'
##' @title Circle distance function
##' @param p A matrix with 2 columns (3 in \code{mesh.dsphere}), each row
##' representing a point in the plane.
##' @param radius radius of circle
##' @param ... additional arguments (not used)
##' @return A vector of length \code{nrow(p)} containing the signed
##' distances to the circle
##' @author Raoul Grasman; translated from original Matlab sources of Per-Olof
##' Persson.
##' @seealso \code{\link{distmesh2d}}, \code{\link{mesh.drectangle}},
##' \code{\link{mesh.diff}}, \code{\link{mesh.intersect}},
##' \code{\link{mesh.union}}
##' @references \url{http://persson.berkeley.edu/distmesh/}
##' 
##' \cite{P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB. SIAM
##' Review, Volume 46 (2), pp. 329-345, June 2004}
##' @keywords arith math
##' @examples
##' 
##' example(distmesh2d)
##' @export
mesh.dcircle <- function (p, radius = 1, ...) {
  if (!is.matrix(p))
    p <- t(as.matrix(p))
  return(sqrt((p^2) %*% c(1, 1)) - radius)
}
