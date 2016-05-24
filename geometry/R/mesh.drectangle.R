##' Signed distance from points \code{p} to boundary of rectangle to
##' allow easy definition of regions in \code{\link{distmesh2d}}.
##' 
##' @title Rectangle distance function
##' @param p A matrix with 2 columns, each row representing a point in
##' the plane.
##' @param x1 lower left corner of rectangle
##' @param y1 lower left corner of rectangle
##' @param x2 upper right corner of rectangle
##' @param y2 upper right corner of rectangle
##' @param \dots additional arguments (not used)
##' @return a vector of length \code{nrow(p)} containing the signed distances
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
##' example(distmesh2d)
##' @export
mesh.drectangle <- function (p, x1 = -1/2, y1 = -1/2, x2 = 1/2, y2 = 1/2, ...)  {
  if (!is.matrix(p)) 
    p = t(as.matrix(p))
  d1 = y1 - p[, 2]
  d2 = -y2 + p[, 2]
  d3 = x1 - p[, 1]
  d4 = -x2 + p[, 1]
  d5 = sqrt(d1^2 + d3^2)
  d6 = sqrt(d1^2 + d4^2)
  d7 = sqrt(d2^2 + d3^2)
  d8 = sqrt(d2^2 + d4^2)
  matmin = function(...) apply(cbind(...), 1, min)
  d = -matmin(matmin(matmin(-d1, -d2), -d3), -d4)
  ix = d1 > 0 & d3 > 0
  d[ix] = d5[ix]
  ix = d1 > 0 & d4 > 0
  d[ix] = d6[ix]
  ix = d2 > 0 & d3 > 0
  d[ix] = d7[ix]
  ix = d2 > 0 & d4 > 0
  d[ix] = d8[ix]
  return(d)
}
