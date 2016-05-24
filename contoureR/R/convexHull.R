#' Get Convex Hull of Points
#' 
#' Returns the sequence of indexes within the supplied numeric vectors \code{x} and \code{y}, that describe the convex 
#' hull containing those points. This is a (slightly modified) implementation of the Andrews Monotone Chain, which is 
#' a well known algorithm that is able  to solve the convex hull with \code{O(nlogn)} complexity. 
#' Typical computation time on a Macbook Air, 1.7Ghz I7, 8Gb Ram, using random points in the range [0,1]:
#' \itemize{ 
#' \item{100K points}{ 0.03 Seconds} 
#' \item{1M points}{ 0.3 seconds, and}
#' \item{10M points}{3.3 seconds.}
#' }
#' @inheritParams getDelaunayMesh
#' @param includeColinear keep or discard the points that lie \strong{ON} the hull, default is to discard (ie dont keep colinear points),
#' as this is the true definition of the convex hull.
#' @rdname getConvexHull
#' @return Returns a vector of integers that represent the '1-based' indexes of the points relative to the 
#' \code{x} and \code{y} input arguments. The resulting vector represents the \strong{closed} list, meaning that the 
#' first index and the last index in the series will be the same.
#' @references https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
#' @examples
#' #Generate the Convex Hull of a Series of Points
#' set.seed(1)
#' x  = runif(100)
#' y  = runif(100)
#' ch = getConvexHull(x,y)
#' 
#' #To demonstrate, Lets view the hull
#' library(ggplot2)
#' df = data.frame(x,y)
#' ggplot(data=df,aes(x,y)) + 
#'    geom_path(data=df[ch,]) + 
#'    geom_point()  
getConvexHull <- function(x,y,includeColinear=FALSE){
  if(!all(is.numeric(x),is.numeric(y))) stop('x and y must be numeric')
  if(length(x) != length(y)) stop('x and y must be the same length')
  convexHullAM_Indexes(x,y,includeColinear=includeColinear,FALSE);
}

