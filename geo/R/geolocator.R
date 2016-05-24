#' Locates points on a plot initialized by geoplot.
#' 
#' The function locates points on a plot initialized by geoplot returning their
#' latitude and longitude.
#' 
#' 
#' @param type Same parameter as in the locator function.  Type = "l" draws
#' line between points.
#' @param n Number of points. Default value is zero, then the point coordinates
#' are located till the mouse' middle button is clicked.
#' @return A list with components \code{$lat} and \code{$lon}.  or (\code{$x,
#' $y} if \code{geopar$projection = "none"})
#' @seealso \code{\link{geoplot}}, \code{\link{locator}},
#' \code{\link{geodefine}}.
#' @export geolocator
geolocator <-
function(type = "p", n = 0)
{
	geopar <- getOption("geopar")
	oldpar <- selectedpar()
	par(geopar$gpar)
	on.exit(par(oldpar))
	if(n == 0)
		x <- locator(type = type)
	else x <- locator(type = type, n = n)
	if(!is.null(x$x)) {
		lat <- invProj(x$x, x$y, geopar$scale, geopar$b0, geopar$b1,
			geopar$l1, projection = geopar$projection)
		if(geopar$projection == "none")
			return(x <- data.frame(x = lat$x, y = lat$y))
		else return(lat <- data.frame(lat = lat$lat, lon = lat$lon))
	}
	else return(list())
}

