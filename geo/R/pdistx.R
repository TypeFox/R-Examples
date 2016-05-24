#' Euclidian distance
#' 
#' Euclidian distance.
#' 
#' 
#' @param y,x Coordinate vectors of starting positions
#' @param y1,x1 Coordinate vectors of end positions
#' @return Euclidian distance
#' @note Check order of coordinates (?), document distance functions together?
#' @seealso Called by \code{\link{variogram}}.
#' @keywords manip
#' @export pdistx
pdistx <-
function(y, x, y1, x1)
{
	return(sqrt((x - x1)^2. + (y - y1)^2.))
}

