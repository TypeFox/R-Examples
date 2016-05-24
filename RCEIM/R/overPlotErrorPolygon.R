overPlotErrorPolygon <- function(x, y, err_y, col="grey", logPlot=FALSE, ...) {	
# A function to superplot an error polygon around the data
	yLower <- (y - err_y)
	if(logPlot) {
		yLower[which(yLower<=0)] <- 1
	}
	yHigher <- (y + err_y)
	xPos <- c(x, x[length(x):1])	
	polygon(xPos, c(yLower, yHigher[length(yHigher):1]), col=col, ...)
}