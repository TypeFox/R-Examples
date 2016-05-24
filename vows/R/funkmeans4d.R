#' Functional k-means clustering for parallel smooths for 4-dimensional data
#' 
#' This is a wrapper function for \code{\link{funkmeans}} to handle 3D image
#' responses.
#' 
#' 
#' @param fdobj a functional data object, of class \code{"\link[fda]{fd}"},
#' defining the set of curves being clustered.
#' @param arr4d a 4-dimensional array containing the raw data that were
#' smoothed at each point.  The first 3 dimensions refer to x, y, and z
#' coordinates and the last dimension corresponds to different images.
#' @param \dots other arguments, passed to \code{\link{funkmeans}}.
#' @return An object of class "funkmeans4d", which is also of class
#' \code{"\link{funkmeans}"} but has the additional component
#' \code{arr.cluster}: an array, of dimension \code{dim(arr4d)[1:3]}, giving
#' the cluster memberships.
#' @author Philip Reiss \email{phil.reiss@@nyumc.org}, Lei Huang
#' \email{huangracer@@gmail.com} and Lan Huo
#' @seealso \code{\link{funkmeans}}, \code{\link{funkpanel}}
#' @examples
#' 
#' # See example for funkpanel
#' @export
funkmeans4d <- function(fdobj, arr4d, ...) {
    has.data = attributes(arr4d)$has.data
    fkmobj <- funkmeans(fdobj=fdobj, ...)
  	# centers.fdobj = fd(coef=fkmobj$fpca$meanfd$coef %*% matrix(1,1,centers) + fkmobj$fpca$harmonics$coef %*% t(fkmobj$centers), basisobj=temp$basis)
  	# include = if (is.null(obj$include)) has.data else obj$include
  	include = has.data
  	arr.cluster = array(NA, dim(include))
  	arr.cluster[include & !is.na(include)] = fkmobj$cluster
  	arr.cluster[has.data & !is.na(include) & !include] = 0
	attr(arr.cluster, "x.ind") = attr(arr4d, "x.ind")
	attr(arr.cluster, "y.ind") = attr(arr4d, "y.ind")
	attr(arr.cluster, "z.ind") = attr(arr4d, "z.ind")
	attr(arr.cluster, "dim.nii") = attr(arr4d, "dim.nii")
	# May need to get some other attributes from arr4d
    # fkmobj$centers.fdobj = centers.fdobj
    fkmobj$arr.cluster = arr.cluster
    fkmobj$include
    class(fkmobj) = c(class(fkmobj), "funkmeans4d")
    fkmobj
}

