
#' Plot a base-graphics scatterplot with accompanying density legend
#' @param x the x coordinates to be handed to plot()
#' @param y the y coordinates of points in the plot()
#' @param colorVar the numeric vector of values used to color the points
#' @param colorPalette a color palette. If 'colorPalette' contains, for example, 
#' 6 colors, then the values of colorVar will be split and assigned to these 6
#' colors
#' @param side the side of the plot to put the density legend on ("left", "right", 
#' "top", or "bottom")
#' @param proportion the proportion of the plot (from 0 to 1) to allocate to the 
#' density legend (defaults to 0.3)
#' @param legendTitle string for labelling the density legend
#' @param ... additional parameters to be passed to plot()
#' @return none, plot is added to device
#' @examples 
#' library(ggplot2)
#' library(RColorBrewer)
#' colorPalette <- brewer.pal(9, "YlOrRd")[4:9]
#' scatterplotDL(x = mtcars$mpg, y = mtcars$wt, colorVar = mtcars$hp, 
#' legendTitle = "Horse Power", colorPalette = colorPalette, pch = 19,
#'		xlab = "MPG (miles per gallon)", ylab = "Weight (tonnes)",  
#'		main = "MPG by Weight in Cars \n Colored by Horse Power")
#' @author Jason Waddell
#' @export
scatterplotDL <- function(x, y, colorVar, colorPalette, 
		side = "right", proportion = 0.3, legendTitle = NULL, ...) {
	
	colorObj <- splitColorVar(colorVar = colorVar, colorPalette)
	colorVec <- colorObj$colorVec
	breaks <- colorObj$breaks
	
	style <- "base"
	
	prepLegend(side = side, proportion = proportion)
	
	if(style == "base") {	
		
		plot(x = x, y = y, col = colorVec, ...)	
		
#	} else if(style == "oa") {
#		
#		oaTemplate(xlim = range(x, na.rm = TRUE), ylim = range(y, na.rm = TRUE), 
#				...)
#		# needs testing... 
#		
	} else {
		stop("Implemented styles are 'oa' and 'base'.")
	}
	
	densityLegend(x = colorVar, colorPalette = colorPalette, side = side, 
			main = legendTitle, colorBreaks = breaks)
}



#' Function to take a numeric vector 'colorVar' and palette 'colorPalette', and 
#' return a list containing a vector of color assignments for each element of
#' 'colorVar' (to be used in plot calls), and a vector of breaks defining the
#' color regions (to be used in densityLegend) 
#' @param colorVar the numeric vector of values used to color the points
#' @param colorPalette a color palette. If 'colorPalette' contains, for example, 
#' 6 colors, then the values of colorVar will be split and assigned to these 6
#' colors
#' @param breaks (optional) a numeric vector of two or more unique cut points 
#' @return  a list containing a vector of color assignments ('colorVec') for each 
#' element of 'colorVar' (to be used in plot calls), and a vector of breaks 
#' ('breaks') defining the color regions (to be used in densityLegend) 
#' @author Jason Waddell
#' @export
splitColorVar <- function(colorVar, colorPalette, breaks = NULL) {
	
	colorVec <- colorPalette[unclass(cut(colorVar, length(colorPalette)))]
	
	if(is.null(breaks))
		breaks <- getBreaks(x = colorVar, breaks = length(colorPalette))
	
	return(list(colorVec = colorVec, breaks = breaks))
}


#' Divide the range of x into intervals, returning the breakpoints of these 
#' intervals
#' @param x a numeric vector which is to be converted to a factor by cutting
#' @param breaks a single number (greater than or equal to 2) giving the number 
#' of intervals into which x is to be cut
#' @param dig.lab integer which is used when labels are not given. It determines
#'  the number of digits used in formatting the break numbers
#' @return a vector of numeric breakpoints
#' @author Jason Waddell
getBreaks <- function (x, breaks, dig.lab = 3L) {
	
	if (!is.numeric(x)) 
		stop("'x' must be numeric")
	if (length(breaks) == 1L) {
		if (is.na(breaks) || breaks < 2L) 
			stop("invalid number of intervals")
		nb <- as.integer(breaks + 1)
		dx <- diff(rx <- range(x, na.rm = TRUE))
		if (dx == 0) {
			dx <- abs(rx[1L])
			breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, 
					length.out = nb)
		} else {
			breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
			breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + 
							dx/1000)
		}
	} else nb <- length(breaks <- sort.int(as.double(breaks)))
	
	if (anyDuplicated(breaks)) 
		stop("'breaks' are not unique")
	
	for (dig in dig.lab:max(12L, dig.lab)) {
		ch.br <- formatC(0 + breaks, digits = dig, width = 1L)
		if (ok <- all(ch.br[-1L] != ch.br[-nb])) 
			break
	}
	
	return(as.numeric(ch.br))
}






#' Create a colored density legend for visually representing the distribution of 
#' a color variable on a plot
#' @param x a numeric vector used to create the density trace
#' @param colorPalette a vector of color values
#' @param colorBreaks a vector of cutoff values for the color regions
#' @param side the side of the plot to place the desntiy legend
#' @param main the main title for the density legend (optional, recommended to 
#' use a title that describes x
#' @return none, graphics are added to the current device
#' @examples 
#' library(ggplot2)
#' library(RColorBrewer)
#' 
#' # subset the data object
#' dsub <- subset(diamonds, x > 5 & x < 6 & y > 5 & y < 6)
#' dsub <- dsub[-which(dsub$z > 4), ]
#' dsub <- dsub[-which(dsub$z < 3), ]
#' 
#' # define color pallette, color vector and color region breaks
#' colorPalette <- brewer.pal(9, "Blues")[4:9]
#' colorObj <- splitColorVar(colorVar = dsub$z, colorPalette)
#' colorVec <- colorObj$colorVec
#' breaks <- colorObj$breaks
#' 
#' # plot the data
#' prepLegend(side = "right", proportion = 0.3)
#' oaTemplate(xlim = range(dsub$x), ylim = range(dsub$y), 
#'		main = "Diamond Length by Width \n Colored by Depth",
#'		xlab = "Length (mm)", ylab = "Width (mm)")
#' points(x = dsub$x, y = dsub$y, col = colorVec, pch = 19, cex = 0.6)
#'
#' # add the legend
#' densityLegend(x = dsub$z, colorPalette = colorPalette, side = "right",
#'		main = "Diamond Depth", colorBreaks = breaks)
#' @author Jason Waddell
#' @export
densityLegend <- function(x, colorPalette, colorBreaks, side = "right",
		main = NULL) {
	
	op <- par('mar')
	par(mar = c(0, 0, 0, 0))
	de1 <- density(x)
	ySpan <- diff(range(de1$x))
	
	if(is.null(colorBreaks)) {
		colorBreaks <- quantile(x, 
				probs = c(0.01, 0.025, 0.16, 0.86, 0.975, 0.99))
		colorBreaks 
	}
	
	if(side %in% c("top", "bottom")){
		ylim <- c(-0.2*max(de1$y), 1.35*max(de1$y))
		xlim <- c(range(de1$x)[1]-0.3*ySpan, range(de1$x)[2]+0.3*ySpan)
	} else {
		xlim <- c(-0.2*max(de1$y), 1.35*max(de1$y))
		ylim <- c(range(de1$x)[1]-0.3*ySpan, range(de1$x)[2]+0.3*ySpan)
	}
	
	plot(x = -1, y = -1, type = "n", ann = FALSE,
			xlim = xlim, ylim = ylim ,	axes = FALSE)
	par(xpd = NA)
	
	
	# plot the bars
	plotBars(de1 = de1, side = side, colorPalette = colorPalette,
			colorBreaks = colorBreaks)
	
	# density trace
	plotDensityTrace(de1 = de1, side = side)
	
	
	# polygon regions
	plotPolygonRegions(de1 = de1, side = side, colorPalette = colorPalette,
			colorBreaks = colorBreaks)
	
	# main
	if(!is.null(main)){
		if(side %in% c("top", "bottom")){
			text(x = xlim[1] + 0.15 * diff(xlim), 
					y = ylim[1] + 0.2 * diff(ylim), main, 
					srt = 90, adj = c(0, 0.5))
		} else {
			text(x = xlim[1] + 0.2 * diff(xlim), 
					y = ylim[1] + 0.85 * diff(ylim), main, adj = c(0, 0.5))
		}
	}
	
	par(mar = op)
}





#' Function to plot all colored density regions of a density legend
#' @param de1 a density() object
#' @param side the side of the plot that the density legend should be plotted on
#' @param colorPalette a vector of color values
#' @param colorBreaks a vector of cutoff values for the color regions
#' @return none, graphics are added to the current device
#' @author Jason Waddell
plotPolygonRegions <- function(de1, side, colorPalette, colorBreaks) {
	
	tempDex <- which(de1$x < colorBreaks[2])
	tempDex <- c(tempDex, max(tempDex) + 1)
	colorPoly(de1, tempDex, col = colorPalette[1], side = side)
	
	for(i in 2:length(colorPalette)){
		tempDex <- which(de1$x > colorBreaks[i] & de1$x < colorBreaks[i+1])
		tempDex <- c(tempDex, max(tempDex) + 1)
		colorPoly(de1, tempDex, colorPalette[i], side = side)
	}
	
	tempDex <- which(de1$x > colorBreaks[i])
	colorPoly(de1, tempDex, col = colorPalette[i], side = side)
	
	if(side %in% c("top", "bottom")){
		lines(x = de1$x, y = 0.3*max(de1$y) + de1$y, col = gray(0.5))
	} else {
		lines(y = de1$x, x = 0.3*max(de1$y) + de1$y, col = gray(0.5))
	}
}


#' Function for plotting a colored polygon as part of a density legend
#' @param de1 a density() object
#' @param tempDex a set of indices corresponding to the range of the current 
#' segment to be plotted (which indices of the density object to ues)
#' @param col the color of the polygon to be plotted
#' @param side the side of the plot that the density legend should be plotted on
#' @return none, graphics are added to the current device
#' @author Jason Waddell
colorPoly <- function(de1, tempDex, col, side){
	
	deRangeValues <- c(de1$x[tempDex], de1$x[max(tempDex)],
			de1$x[min(tempDex)], de1$x[min(tempDex)])
	deDepthValues <- 0.3*max(de1$y) + 
			c(de1$y[tempDex], 0, 0, de1$y[min(tempDex)])
	
	if(side %in% c("top", "bottom")){
		polygon(x = deRangeValues, y = deDepthValues, col = col, border = NA)
	} else {
		polygon(y = deRangeValues, x = deDepthValues, col = col, border = NA)
	}
}




#' Function for plotting the density trace outline in a density legend
#' @param de1 a density() object
#' @param side the side of the plot that the density legend should be plotted on
#' @return none, graphics are added to the current device
#' @author Jason Waddell
plotDensityTrace <- function(de1, side) {
	
	if(side %in% c("top", "bottom")){
		segments(y0 = 0.3*max(de1$y), x0 = min(de1$x), 
				x1 = max(de1$x), col = gray(0.8))	
		lines(x = de1$x, y = 0.3*max(de1$y) + de1$y, col = gray(0.5))
	} else {
		segments(x0 = 0.3*max(de1$y), y0 = min(de1$x), 
				y1 = max(de1$x), col = gray(0.8))	
		lines(y = de1$x, x = 0.3*max(de1$y) + de1$y, col = gray(0.5))
	}	
	
}





#' A function for creating the segmented color bars in a density legend
#' @param de1 a density() object
#' @param side the side of the plot that the density legend should be plotted on
#' @param colorPalette A vector of color values
#' @param colorBreaks A vector of cutoff values for the color regions
#' @return none, graphics are added to the current device
#' @author Jason Waddell
plotBars <- function(de1, side, colorPalette, colorBreaks) {
	
	yloc <- pretty(de1$x)
	yloc <- yloc[yloc <= max(de1$x) & yloc >= min(de1$x) ]
	
	if(side %in% c("top", "bottom")){
		
		# Axis code		
		xloc <- c(-0.1, 0.05)*max(de1$y)
		
		segments(y0 = xloc[1], y1 = xloc[2], x0 = yloc, col = gray(0.6), lwd = 1.5)
		text(x = yloc - 0.006*diff(range(yloc)), y = xloc[1], labels = yloc, 
				adj = c(1, 0), col = gray(0.5)) 
		
		# bars
		xleft = 0.1*max(de1$y); xright = 0.2*max(de1$y)
		rect(ybottom = xleft, ytop = xright, 
				xleft = min(de1$x), xright = max(de1$x), col = colorPalette[1], 
				border = FALSE)
		
		for(i in 2:(length(colorPalette)-1))
			rect(ybottom = xleft, ytop = xright, 
					xleft = colorBreaks[i], xright = colorBreaks[i+1], col = colorPalette[i], 
					border = FALSE)		
		
		rect(ybottom = xleft, ytop = xright, 
				xleft = colorBreaks[i+1], xright = max(de1$x), 
				col = colorPalette[i+1], border = FALSE)
		
	} else {   # left, right
		
		# Axis code		
		xloc <- c(-0.1, 0.05)*max(de1$y)
		
		# tick marks
		segments(x0 = xloc[1], x1 = xloc[2], y0 = yloc, col = gray(0.6), lwd = 1.5)
		text(y = yloc + 0.006*diff(range(yloc)), x = xloc[2], labels = yloc, 
				adj = c(1, 0), col = gray(0.5)) 
		
		# bars
		xleft = 0.1*max(de1$y); xright = 0.2*max(de1$y)
		rect(xleft = xleft, xright = xright, 
				ybottom = min(de1$x), ytop = max(de1$x), col = colorPalette[1], 
				border = FALSE)
		
		for(i in 2:(length(colorPalette)-1))
			rect(xleft = xleft, xright = xright, 
					ybottom = colorBreaks[i], ytop = colorBreaks[i+1], col = colorPalette[i], 
					border = FALSE)		
		
		rect(xleft = xleft, xright = xright, 
				ybottom = colorBreaks[i+1], ytop = max(de1$x), 
				col = colorPalette[i+1], border = FALSE)
	}
}









##' Create a density legend plot
##' @param x vector of values with which to create a density
##' @param breaks vector of breakpoints for the density, from the quantile() function
##' @param colVec a vector of colors, of length length(breaks)-1
##' @param main (optional) a title for the density legend
##' @return none. plot is added to the current device
##' @author Jason Waddell
#densityLegendDecreped <- function(x, breaks, colVec, main){
#	de1 <- density(x)
#	ySpan <- diff(range(x, na.rm = TRUE))
#	x <- de1$x; y <- de1$y
#	
#	op <- par(no.readonly = TRUE)
#	par(plt = c(0.05, 0.95, 0.15, 0.85))
#	par(xpd = NA)
#	blankPlot( xlim = c(-0.2*max(y),1.35*max(y)), 
#			ylim = c(range(x)[1]-0.25*ySpan, range(x)[2]+0.25*ySpan)  )
#	
#	xleft = 0; xright = 0.2*max(y)
#	for(i in 1:(length(breaks)-1))
#		rect(xleft = xleft, xright = xright, 
#				ybottom = breaks[i], ytop = breaks[i+1], col = colVec[i], border = FALSE)		
#	rect(xleft = xleft, xright = xright, 
#			ybottom = min(x), ytop = breaks[1], col = colVec[1], border = FALSE)
#	rect(xleft = xleft, xright = xright, 
#			ybottom = max(x), ytop = breaks[length(breaks)], col = colVec[length(colVec)], border = FALSE)
#	
#	if(!missing(main))
#		text(x = xleft, y = max(x) + 0.03*ySpan, labels = main, cex = 1.2, col = gray(0.5), adj = c(0, 0))
#	
#	segments(x0 = 0.3*max(y), y0 = min(x), 
#			y1 = max(x), col = gray(0.8))	
#	lines(y = x, x = 0.3*max(y) + y, col = gray(0.5))
#	
#	
#	colorPoly <- function(tempDex, col, minX, maxX){
#		if(missing(minX) & !missing(maxX)){
#			propTop <- (maxX - x[max(tempDex)])/(x[max(tempDex)+1] - x[max(tempDex)])
#			yTop <- y[max(tempDex)] + propTop*(y[max(tempDex)+1] - y[max(tempDex)])
#			
#			polygon(y = c(x[tempDex], maxX, maxX, x[min(tempDex)], x[min(tempDex)]),
#					x = 0.3*max(y) + c(y[tempDex], yTop, 0, 0, y[min(tempDex)]), col = col, 
#					border = NA	)
#		}
#		
#		if(!missing(minX) & !missing(maxX)){
#			propTop <- (maxX - x[max(tempDex)])/(x[max(tempDex)+1] - x[max(tempDex)])
#			yTop <- y[max(tempDex)] + propTop*(y[max(tempDex)+1] - y[max(tempDex)])
#			
#			propBottom <- (x[min(tempDex)] - minX)/(x[min(tempDex)] - x[min(tempDex)-1])
#			yBottom <- y[min(tempDex)] - propBottom*(y[min(tempDex)] - y[min(tempDex) - 1])
#			
#			polygon(y = c(x[tempDex], maxX, maxX, minX, minX),
#					x = 0.3*max(y) + c(y[tempDex], yTop, 0, 0, yBottom), col = col, 
#					border = NA	)
#		}
#		
#		if(!missing(minX) & missing(maxX)){
#			propBottom <- (x[min(tempDex)] - minX)/(x[min(tempDex)] - x[min(tempDex)-1])
#			yBottom <- y[min(tempDex)] - propBottom*(y[min(tempDex)] - y[min(tempDex) - 1])
#			
#			polygon(y = c(x[tempDex], x[max(tempDex)], x[max(tempDex)], minX, minX),
#					x = 0.3*max(y) + c(y[tempDex], y[max(tempDex)], 0, 0, yBottom), col = col, 
#					border = NA	)
#		}
#	}
#	
#	tempDex <- which(x <= breaks[2])
#	colorPoly(tempDex, col = colVec[1], maxX = breaks[2])
#	
#	for(i in 2:(length(breaks)-1)){
#		tempDex <- which(x <= breaks[i+1] & x >= breaks[i])
#		colorPoly(tempDex, colVec[i], minX = breaks[i], maxX = breaks[i+1])
#	}
#	
#	tempDex <- which(x > breaks[length(breaks)-1])
#	colorPoly(tempDex, col = colVec[length(colVec)], minX = breaks[length(breaks)-1])
#	
#	lines(y = x, x = 0.3*max(y) + y, col = gray(0.5))
#	
#	
#	# Axis code
#	par(xpd = NA)
#	yloc <- pretty(x)
#	yloc <- yloc[yloc <= max(x) & yloc >= min(x) ]
#	
#	
#	xloc <- c(-0.25, -0.05)*max(y)
#	
#	segments(x0 = xloc[1], x1 = xloc[2], y0 = yloc, col = gray(0.8), lwd = 1.5)
#	text(y = yloc + 0.006*diff(range(yloc)), x = xloc[2], labels = yloc, adj = c(1, 0), col = gray(0.8))
#	par(op)
#}
#
