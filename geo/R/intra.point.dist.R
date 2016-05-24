#' Intra point/position distance
#' 
#' For a data frame of positions return the vector of intra point distances
#' 
#' 
#' @param x List of positions with components \code{lat} and \code{lon}.
#' @return Vector of distances between the points in \code{x}.
#' @seealso \code{\link{arcdist}}
#' @keywords arith
#' @examples
#' 
#' # distances along the perimeter of a statistical rectangle
#' pos <- rPeri(323)
#' intra.point.dist(pos)
#' sum(intra.point.dist(pos))
#' 
#' @export intra.point.dist
intra.point.dist <-
function(x)
{
	n <- length(x$lat)
	arcdist(x$lat[ - n], x$lon[ - n], x$lat[-1], x$lon[-1])
}

