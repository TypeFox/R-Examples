#' Add arrowhead to plot
#' 
#' Adds arrowhead at the start or the end of a list of coordinates on a
#' geoplot.
#' 
#' 
#' @param data Dataframe with coordinates in columns \code{lat, lon}.
#' @param start Start from beginning (i.e. arrow points at first coordinate),
#' default TRUE.
#' @param size Size of arrow
#' @param \dots Additional arguments to \code{Arrow}
#' @return Draws an Arrowhead at the start or the end of the coordinates given
#' in \code{data}, returns invisibly the outline as four coordinate pairs in a
#' dataframe.
#' @note May need further elaboration and/or detail.
#' @seealso Called by \code{\link{geocurve}}, calls \code{\link{findline}},
#' \code{\link{invProj}} and \code{\link{Arrow}}.
#' @keywords aplot
#' @export geolines.with.arrows
geolines.with.arrows <-
function(data, start = T, size = 0.2, ...)
{
	geopar <- getOption("geopar")
	if(!is.data.frame(data))
		data <- data.frame(data)
	n <- nrow(data)
	if(start)
		i <- c(1:n)
	else i <- seq(n, 1, by = -1)
	tmpdata <- data[i,  ]
	limits <- invProj(geopar$limx, geopar$limy)
	plt.size <- geopar$gpar$pin
	dlon <- size/plt.size[1] * diff(limits$lon)
	dlat <- size/plt.size[2] * diff(limits$lat)
	theta <- seq(0, 2 * pi, by = 0.1)
	lat <- tmpdata[1, "lat"] + dlat * sin(theta)
	lon <- tmpdata[1, "lon"] + dlon * cos(theta)
	circle <- data.frame(lat = lat, lon = lon)
	xr <- findline(tmpdata, circle, plot = F)
	i <- is.na(xr$lat)
	i1 <- c(1:length(i))
	i1 <- i1[i]
	n <- min(i1) - 1
	pos <- list(lat = c(xr$lat[n], tmpdata$lat[1]), lon = c(xr$lon[n],
		tmpdata$lon[1]))
	pos <- Arrow(pos, ...)
	return(invisible(pos))
}

