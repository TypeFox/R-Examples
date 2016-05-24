#' Find a circle around a point
#' 
#' Through simple approximate geometry the perimeter of a circle around a point
#' given in lat and lon is found
#' 
#' 
#' @param lat,lon Latitude and longitude in decimal degrees.
#' @param rad radius of the circle in \code{nmi}
#' @param n number of points on the circle perimeter
#' @keywords hplot
#' @examples
#' 
#' ## draws a circle around Iceland
#' geoplot()
#' geolines(circle.one.point(65,-19,150,n=50),lwd=3,col="red")
#' 
#' @export circle.one.point
circle.one.point <-
function(lat, lon = NULL, rad, n = 10.)
{
	if(is.null(lon)) {
		lon <- lat$lon
		lat <- lat$lat
	}
	out <- list(lat = numeric(n), lon = numeric(n))
	dlat <- rad/60.
	dlon <- rad/arcdist(lat, lon, lat, lon - 1.)
	angles <- seq( - pi, pi, length = n)
	out$lat <- lat + dlat * sin(angles)
	out$lon <- lon + dlon * cos(angles)
	data.frame(out)
}

