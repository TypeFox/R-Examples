#' Adds points on plots initialized by geoplot.
#' 
#' Plot points on a graph initialized by geoplot. Data is stored as lat, lon or
#' x,y depending on the projection. The program plots the transformation of the
#' data. Parameters for the projection are stored in the list geopar. Similar
#' to the Splus function points.
#' 
#' 
#' @param lat,lon Plot points on a graph initialized by geoplot. Data is stored
#' as lat, lon or x,y depending on the projection. The program plots the
#' transformation of the data. Parameters for the projection are stored in the
#' list geopar. Similar to the Splus function points.
#' @param pch Type of symbol used options are for example. " ","*","+","." or
#' anything else, letter, digit, or symbol. Default is "*"
#' @param cex Relative size of character and symbols (see the help on the
#' parameter cex).  The size of plotted characeters is cex time the parameter
#' csi that can be seen by \code{par("csi")}.  In earlier versions of geoplot
#' the parameter csi was set but csi is a parameter that can not be set in R.
#' The parameter mkh should probably be used for symbols instead of cex, see
#' help on graphical parameters.
#' @param col Colour number used.  Default value is one.
#' @param lwd Linewidth used. Default is the value set when the program was
#' called.
#' @param outside If TRUE geopoints will plot points outside the specified
#' limits set by geoplot(). If FALSE, which is default, outside points will be
#' skipped.
#' @param jitter useful if many datapoints have the same coordinates, points
#' are jittered randomly to make common values look bigger.jitter=0.016 is
#' often about right but you may want to have jitter smaller or bigger varying
#' on plot.
#' @param mkh Size of symbol in inches.  If not given cex is used instead.
#' @param csi Size of character.  This parameter can not be set in R but for
#' compatibility with old Splus scripts the parameter cex is readjusted by
#' \code{cex = cex*csi/0.12}.  Use of this parameter is not recommended.
#' Default value is NULL i.e not used.
#' @seealso \code{\link{geoplot}}, \code{\link{geopolygon}},
#' \code{\link{geolines}}, \code{\link{points}}, \code{\link{geotext}},
#' \code{\link{geosymbols}}, \code{\link{geocontour.fill}},
#' \code{\link{geogrid}}, \code{\link{geocontour}}.
#' @examples
#' 
#' \dontrun{       geopoints(deg)                  # Plots * in the points
#'                                                 # defined by deg$lat,deg$lon.
#' 
#'        geopoints(deg$lat,deg$lon,pch="*",col=5) # Same but uses color 5.
#' 
#'        geopoints(fd$x,fd$y)                     # Points in x,y when
#'                                                 # projection in geoplot
#'                                                 # was "none".
#' }
#' @export geopoints
geopoints <-
function(lat, lon = 0, pch = "*",cex =0.7, col = 1, lwd = 0, outside = FALSE,
	jitter = NULL, mkh = NULL,csi=NULL)
{
   geopar <- getOption("geopar")
   if(!is.null(csi)) cex <- cex*csi/0.12  # Compatibility with old program
	if(length(lon) == 1 && length(lat) > 1) {
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
	if(!is.null(jitter)) {
		lat <- lat + runif(length(lat), -1, 1) * jitter
		lon <- lon + runif(length(lon), -1, 1) * jitter * 2
	}
	oldpar <- selectedpar()
	par(geopar$gpar)
	if(lwd != 0)
		par(lwd = lwd)
	if(outside)
		par(xpd = TRUE)
	else par(xpd = FALSE)
	on.exit(par(oldpar))
	par(cex = cex)
	xx <- Proj(lat, lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
		geopar$projection)
	if(!outside) {
		ind <- c(1:length(lat))
		ind <- ind[xx$x > geopar$limx[2] | xx$x < geopar$limx[1] | xx$
			y < geopar$limy[1] | xx$y > geopar$limy[2] | is.na(
			xx$x) | is.na(xx$y)]
		if(length(ind) > 0) {
			xx$x <- xx$x[ - ind]
			xx$y <- xx$y[ - ind]
		}
	}
	if(!is.null(mkh))
		points(xx$x, xx$y, pch = pch, col = col, mkh = mkh)
	else points(xx$x, xx$y, pch = pch, col = col)
	return(invisible())
}

