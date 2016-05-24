#' Identifies points on plots using lat and lon coordinates.
#' 
#' Works the same way as identify except that it also accepts coordinates as
#' lat and lon.  Identifies points on a plot identified by the user.
#' 
#' Observations that have missing values in either lat or lon are treated as if
#' they were not given.  When using the X11 driver under the X Window System, a
#' point is identified by positioning the cursor over the point and pressing
#' the left button.  To exit identify press the middle button (both buttons on
#' a two button mouse) while the cursor is in the graphics window.  The same
#' procedure is used under the suntools driver.  This function may also be used
#' with the "tek" drivers.
#' 
#' Some devices that do not allow interaction prompt you for an x,y pair. The
#' nearest point to the locator position is identified, but must be at most 0.5
#' inches away.  In case of ties, the earliest point is identified.
#' 
#' @param lat,lon Coordinates of points.  The coordinates can be given by two
#' vectors or a data.frame containin vectors \code{lat} and \code{lon}.
#' @param labels Vector giving labels for each of the points.  If supplied,
#' this must have the same length as lat and lon.  As a default the vector
#' indece number of the points will be used.
#' @param n Maximum number of points to be identified.
#' @param plot If true, geoidentify plots the labels if the points identified.
#' @param atpen If true, plotted identification is relative to locator position
#' when the point is identified; otherwise, plotting is relative to the
#' identified lat,lon value. This can be useful when points are crowded.
#' Default is true.
#' @param offset Identification is plotted as a text string, moved offset
#' charecter from the point.  If the locator was left (right) of the nearest
#' point, the label will be offset to the left (right) of the point.
#' @param col The color of the labels.
#' @param cex Character size expansion of label characters.
#' @return Indeces (in lat and lon) corresponding to the identified points.
#' @section Side Effects: Labels are placed on the current plot if plot is
#' true.
#' @seealso \code{\link{identify}}, \code{\link{geolocator}},
#' \code{\link{geotext}}.
#' @examples
#' 
#' \dontrun{       geoidentify(stations, labels = stations$temp)
#'        # plots the temperature in the closest measuring point.
#' 
#'        geoidentify(stations, atpen = FALSE) 
#'        # plots the indece number of the station closest to
#'        # where pointed at the stations position.
#' }
#' @export geoidentify
geoidentify <-
function(lat, lon = NULL, labels = 1, n = 0, plot = TRUE, atpen = TRUE, offset = 0.5,
	col = 1, cex = 1)
{
	geopar <- getOption("geopar")
	oldpar <- selectedpar()
	par(geopar$gpar)
	par(cex = cex)
	par(col = col)
	on.exit(par(oldpar))
	if(is.null(lon)) {
		if(geopar$projection == "none") {
			lon <- lat$y
			lat <- lat$x
		}
		else {
			lon <- lat$lon
			lat <- lat$lat
		}
	}
	if(geopar$projection != "none") {
		# degrees and minutes
		if(mean(lat, na.rm = TRUE) > 1000) {
			lat <- geoconvert(lat)
			lon <-  - geoconvert(lon)
		}
	}
	if(length(labels) == 1 && length(lat) > 1)
		labels <- seq(along = lat)
	if(n == 0)
		n <- length(lat)
	xx <- Proj(lat, lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
		geopar$projection)
	z <- identify(xx$x, xx$y, labels = labels, n = n, plot = plot, atpen = 
		atpen, offset = offset)
	return(z)
}

