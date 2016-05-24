##' Signed distance from points \code{p} to boundary of sphere to
##' allow easy definition of regions in \code{\link{distmeshnd}}.
##'
##' @title Sphere distance function
##' @param p A matrix with 2 columns (3 in \code{mesh.dsphere}), each row
##' representing a point in the plane.
##' @param radius radius of sphere
##' @param ... additional arguments (not used)
##' @return A vector of length \code{nrow(p)} containing the signed
##' distances to the sphere
##' @author Raoul Grasman; translated from original Matlab sources of Per-Olof
##' Persson.
##' @seealso \code{\link{distmeshnd}}
##' @references \url{http://persson.berkeley.edu/distmesh/}
##' 
##' \cite{P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB. SIAM
##' Review, Volume 46 (2), pp. 329-345, June 2004}
##' @keywords arith math
##' @examples
##' 
##' example(distmeshnd)
##' @export
"mesh.dsphere" <-
function (p, radius = 1, ...) 
{
    if (!is.matrix(p)) 
        p = t(as.matrix(p))
    sqrt((p^2) %*% rep(1, ncol(p))) - radius
}
