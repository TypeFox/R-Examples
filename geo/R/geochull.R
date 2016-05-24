#' Convex hull of a set of positions
#' 
#' Finds the convex hull of a set of positions in lat and lon.
#' 
#' 
#' @param lat,lon Position(s) as decimal degrees latitude and longitude.  If
#' 'lat' is 'list' its components 'lat$lat' and 'lat$lon' are used for 'lat'
#' and 'lon'.
#' @return List with components \code{lat} and \code{lon} of the convex hull of
#' the input.
#' @seealso \code{\link{Proj}}, \code{\link{invProj}},
#' \code{\link[grDevices]{chull}}
#' @keywords manip
#' @examples
#' 
#' # draws the convex hull of Iceland's coastline.
#' geoplot()
#' geolines(geochull(island))
#' 
#' @export geochull
geochull <-
function(lat, lon = NULL)
{
	if(is.null(lon)) {
		lon <- lat$lon
		lat <- lat$lat
	}
	x <- Proj(lat, lon)$x
	y <- Proj(lat, lon)$y
	id <- chull(x, y)
	id <- c(id, id[1])
	data.frame(invProj(x[id], y[id])[c("lat", "lon")])
}

