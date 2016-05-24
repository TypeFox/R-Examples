#' Smooth curves and put arrows in the beginning and end.
#' 
#' Smooth curves and put arrows in the beginning and end.  Useful for plotting
#' closed areas, making smooth lines with arrows etc.  The function needs
#' dataframe with columns lat and lon.  No NA values currently allowed. Uses
#' the functions ns or ps to smooth the curve depending on if the curve is open
#' or closed.  Smooths seperatly lat ~ ns(d,df) , lon ~ ns(d,df).  d is here
#' distance along the curve.  If the same curve is smoothed many times the
#' output is stored the first time and the function called later with smooth=F.
#' 
#' 
#' @param data Dataframe that must have columns lat and lon.
#' @param df Degrees of Freedom in smoothing the curve.  Default is number of
#' points/2.  Maximum allowed df is nrow(data)-1.  df="all" gives this value.
#' @param n How much denser the output data is than the input data.  n=10 means
#' that the value of the smoothing spline is predicted at 9 evenly spaced
#' points beetween each pair of data points.  Default is n = 10.
#' @param open Is the curve open or closed.  Default open.
#' @param arrow If there is an arrow in the beginning, end or both of the line
#' segment.  Options are "none","start","end","both".
#' @param col Color of line.
#' @param lwd Width of line segment. Default is 1.  lwd=2 seems to fit well
#' with default size of arrows.
#' @param size Size of arrow, default 0.2 inches.
#' @param angle Angle of arrow opening, default 15 degrees.
#' @param smooth Should the data be smoothed, default is T.
#' @param plot Should the curve be plotted, default is T.
#' @param \dots Other parameters to geolines and geopolygons.
#' @return If smooth=F = the input data is returned, else smoothed inputdata.
#' @section Side Effects: Plots the line segment on the screen.
#' @seealso \code{\link{geolines}}, \code{\link{geopolygon}}, \code{\link{ns}}.
#' @examples
#' 
#' \dontrun{       # define curve store the result.
#'        curve1 <- geocurve(geolocator(type="p"),arrow="end",lwd=2)
#'        # use the result.
#'        geocurve(curve1,smooth=F,arrow="start",lwd=2,col=150)
#' 
#'        # define closed area and hatch it.
#'        area1 <- geocurve(geodefine(),open=F)
#'        geopolygon(area1,density=10,col=1)
#' 
#'        # Make closed curve with big arrow and not store the result.
#'        geocurve(geodefine(),open=F,arrow="end",lwd=2,size=0.5)
#' }
#' @export geocurve
geocurve <-
function(data, df = nrow(data)/2, n = 10, open = T, arrow = "none", col = 1,
	lwd = 1, size = 0.2, angle = 15, smooth = T, plot = T, ...)
{
	if(df == "all")
		df <- nrow(data) - 1
	if(smooth) {
		if(open)
			data <- Open.curve(data, df, n)
		else {
			n <- nrow(data)
			if(data[1, "lat"] != data[n, "lat"] || data[1, "lon"] !=
				data[n, "lon"])
				data <- rbind(data, data[1,  ])
			data <- Closed.curve(data, df, n)
		}
	}
	if(plot) {
		geolines(data, lwd = lwd, col = col, ...)
		if(arrow == "start" || arrow == "both")
			geolines.with.arrows(data, start = T, col = col, size
				 = size, angle = angle)
		if(arrow == "end" || arrow == "both")
			geolines.with.arrows(data, start = F, col = col, size
				 = size, angle = angle)
	}
	return(invisible(data))
}

