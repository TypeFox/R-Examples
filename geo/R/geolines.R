#' Add lines to current plot initialized by geoplot.
#' 
#' Add lines to a plot initialized by geoplot. Data is stored as lat, lon or
#' x,y. Lists are assumed to have the components \code{$x} and \code{$y} if
#' projection in geoplot was "none", else \code{$lat},\code{$lon}. The program
#' transforms the data as specified in geoplot. Similar to the Splus function
#' lies.
#' 
#' 
#' @param lat Latitude of data. ( or x coordinate)
#' @param lon Longitude of data. ( or y coordinate) Negative values mean
#' western longitudes. Default value is zero. If lon is zero then the data is
#' stored as \code{lat$lat} and \code{lat$lon}. (or \code{lat$x} and
#' \code{lat$y})
#' @param col Colour number used for plotting the lines, default value is 1.
#' @param lwd Line width. Default is to use the width set when the program was
#' called.
#' @param lty Line type. Default is to use the width set when the program was
#' called. See Splus manuals for numbers corresponding to different linetypes
#' and linewidths.
#' @param nx Parameter only used with Lambert transform when lines in lat,lon
#' are curves in x,y. If nx > 1, nx-1 points are put between each two
#' datapoints in lat, lon before projection is done. For example:
#' 
#' \code{geolines(c(66, 66), c(-30, -10), nx = 50)}
#' 
#' plots a line onto the 66 degree latitude from -30 to -10. The line is curved
#' because it is made of 50 segments.
#' @param outside Logical, should lines outside the plot region be drawn?
#' Default FALSE.
#' @param return.data Logical, should the data be returned? Default FALSE.
#' @return No values returned.
#' @section Side Effects: The projection is stored in geopar\$projection and
#' parameters for the transform in \code{geopar$b0, geopar$b1 and geopar$l1}.
#' @seealso \code{\link{geoplot}}, \code{\link{geopolygon}},
#' \code{\link{geopoints}}, \code{\link{geotext}}, \code{\link{geosymbols}},
#' \code{\link{geocontour.fill}}, \code{\link{geogrid}},
#' \code{\link{geocontour}}.
#' @keywords aplot
#' @examples
#' 
#'        geolines(island)                      # plot iceland.
#'        geolines(island$lat, island$lon, col = 1) # same.
#' 
#'        #######################################################
#' 
#'        geoplot(xlim=c(0, -50), ylim=c(60, 75), projection = "Lambert")
#'        # Set up a Lambert plot.
#' 
#'        geolines(c(66, 66), c(-30, -10), nx = 50, col = 155, lwd = 2)
#'        # Draw a line with colour 155 and width 2.
#' 
#'        geopolygon(island)
#'        geolines(island, col = 3, lwd = 3)
#'        geolines(eyjar, col = 40)
#'        geolines(faeroes, col = 40)
#'        geolines(greenland, col = 3, lwd = 3)
#' #       geolines(janmayen, col = 40)
#'        # Plot some more countries using geolines.
#' 
#' @export geolines
geolines <-
function(lat, lon = 0, col = 1, lwd = 0, lty = 0, nx = 1, outside = FALSE, 
	return.data = FALSE)
{
	geopar <- getOption("geopar")
	if(length(lon) == 1) {
		# For polygon structures.
		if(!is.null(lat$length)) n <- lat$length else n <- max(c(length(
				lat$y), length(lat$lat)))
		if(geopar$projection == "none") {
			lon <- lat$y[1:n]
			lat <- lat$x[1:n]
		}
		else {
			lon <- lat$lon[1:n]
			lat <- lat$lat[1:n]
		}
	}
	if(geopar$projection != "none") {
		# degrees and minutes
		if(mean(lat, na.rm = TRUE) > 1000) {
			lat <- geoconvert(lat)
			lon <-  - geoconvert(lon)
		}
	}
	if(outside)
		par(xpd = TRUE)
	else par(xpd = FALSE)
	if(nx > 1) {
		# fill in with points for lambert. 
		x <- fill.points(lat, lon, nx, option = 2)
		lat <- x$x
		lon <- x$y
	}
	oldpar <- selectedpar()
	par(geopar$gpar)
	if(lwd != 0)
		par(lwd = lwd)
	if(lty != 0)
		par(lty = lty)
	on.exit(par(oldpar))
	xx <- Proj(lat, lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
		geopar$projection)
	if(!outside) {
		gx <- geopar$limx
		gy <- geopar$limy
		border <- list(x = c(gx[1], gx[2], gx[2], gx[1], gx[1]), y = c(
			gy[1], gy[1], gy[2], gy[2], gy[1]))
		xx <- findline(xx, border)
	}
	else par(xpd = FALSE)
	#c program not used
	lines(xx$x, xx$y, col = col)
	par(oldpar)
	if(return.data) {
		xx <- invProj(xx)
		xx <- data.frame(lat = xx$lat, lon = xx$lon)
		return(invisible(xx))
	}
	else return(invisible())
}

