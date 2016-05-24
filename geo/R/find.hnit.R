#' Find coordinate(s) ??
#' 
#' Find coordinates with some sort of interpolation (??).
#' 
#' 
#' @param pt Point(s) ??
#' @param poly Polygon (??)
#' @return Returns list with components: \item{x, y }{of coordinates}
#' @note ~~further notes~~ Needs elaboration.
#' @seealso Called by \code{\link{findcut}}
#' @keywords manip
#' @export find.hnit
find.hnit <-
function(pt, poly)
{
	pt1 <- floor(pt)
	pt2 <- pt - pt1
	y <- poly$y[pt1] + pt2 * (poly$y[pt1 + 1] - poly$y[pt1])
	x <- poly$x[pt1] + pt2 * (poly$x[pt1 + 1] - poly$x[pt1])
	return(data.frame(x = x, y = y))
}

