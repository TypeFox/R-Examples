##' Compute the signed distances from points \code{p} to a region
##' defined by the difference, union or intersection of regions
##' specified by the functions \code{regionA} and \code{regionB}.
##' \code{regionA} and \code{regionB} must accept a matrix \code{p}
##' with 2 columns as their first argument, and must return a vector
##' of length \code{nrow(p)} containing the signed distances of the
##' supplied points in \code{p} to their respective regions.
##'
##' @title Difference, union and intesection operation on  two regions
##' @return A vector of length \code{nrow(p)} containing the signed
##' distances to the boundary of the region.
##' @author Raoul Grasman; translated from original Matlab sources of Per-Olof
##' Persson.
##' @param p A matrix with 2 columns (3 in \code{mesh.dsphere}), each row
##' representing a point in the plane.
##' @param regionA vectorized function describing region A in the
##' union / intersection / difference
##' @param regionB vectorized function describing region B in the
##' union / intersection / difference
##' @param ... additional arguments passed to \code{regionA} and
##' \code{regionB}
##' @aliases mesh.diff mesh.union mesh.intersect
##' @seealso \code{\link{distmesh2d}}, \code{\link{mesh.dcircle}},
##' \code{\link{mesh.drectangle}} \code{\link{mesh.dsphere}}
##' @export mesh.diff mesh.union mesh.intersect
mesh.diff <-function (p, regionA, regionB, ...) {
  return(matmax(regionA(p, ...), -regionB(p, ...)))
}
