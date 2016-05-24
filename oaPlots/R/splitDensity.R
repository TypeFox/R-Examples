#' Utility Function to Draw a Split Density 
#' @param x x vector from a density object. e.g. data <- rnorm(100); x <- density(data)$x
#' @param y y vector from a density object
#' @param yshift vertical shift to be applied to the y object
#' @param colVec shade for the interior of the polygon
#' @param colOut the color for the outer density line
#' @param lwd line width for the outer density line
#' @param split vector specifying color splits
#' @param i the panel number
#' @return none. polygon is added to the current device
#' @author Jason Waddell
#' @noRd
drawPolyMiddle <- function(x, y, yshift = 0, colShade, colOut, lwd = 2, split, i){
	polygon(c(split[i-1], split[i-1], x, split[i], split[i], split[i-1], split[i-1]), 
			c(0, y[1], y, y[length(y)], 0, 0, y[1])+yshift, border = NA, col = colShade)
}


#' Utility Function to Draw a Split Density
#' @param x x vector from a density object. e.g. data <- rnorm(100); x <- density(data)$x
#' @param y y vector from a density object
#' @param yshift vertical shift to be applied to the y object
#' @param colVec shade for the interior of the polygon
#' @param colOut the color for the outer density line
#' @param lwd line width for the outer density line
#' @param split value specifying where to place the leftmost color split
#' @return none. polygon is added to the current device
#' @author Jason Waddell
#' @noRd
drawPolyS1 <- function(x, y, yshift = 0, colShade, colOut, lwd = 2, split){
	polygon(c(x, split, split, x[1], x[1]), 
			c(y, y[length(y)], 0, 0, y[1])+yshift, border = NA, col = colShade)
}

#' Utility Function to Draw a Split Density 
#' @param x x vector from a density object. e.g. data <- rnorm(100); x <- density(data)$x
#' @param y y vector from a density object
#' @param yshift vertical shift to be applied to the y object
#' @param colVec shade for the interior of the polygon
#' @param colOut the color for the outer density line
#' @param lwd line width for the outer density line
#' @param split value specifying where to place the rightmost color split
#' @return none. polygon is added to the current device
#' @author Jason Waddell
#' @noRd
drawPolyS2 <- function(x, y, yshift = 0, colShade, colOut, lwd = 2, split){
	polygon(c(x, x[length(x)], split, split), 
			c(y, 0, 0, y[1])+yshift, border = NA, col = colShade)
}

#' Utility Function to Draw a Split Density
#' @param x x vector from a density object. e.g. data <- rnorm(100); x <- density(data)$x
#' @param y y vector from a density object
#' @param yshift vertical shift to be applied to the y object
#' @param colVec shade for the interior of the polygon
#' @param colOut the color for the outer density line
#' @param lwd line width for the outer density line
#' @return none. polygon is added to the current device
#' @author Jason Waddell
#' @noRd
drawPoly <- function(x, y, yshift = 0, colShade, colOut, lwd = 2){
	polygon(c(x, x[length(x)], x[1], x[1]), c(y, 0, 0, y[1])+yshift, border = NA, col = colShade)
	lines(x, y+yshift, col = colOut, lwd = lwd)
}


#' Draw a Split Density Plot
#' @param x x vector from a density object. e.g. data <- rnorm(100); x <- density(data)$x
#' @param y y vector from a density object
#' @param densityObj an object created by the function density()
#' @param yshift vertical shift to be applied to the y object
#' @param colVec color vector for the shaded regions that compose the interior of the plot. 
#' The length of 'colVec' should be one greater than the length of split
#' @param outerCol the color for the outer density line
#' @param lwd line width for the outer density line
#' @param split vector of x values at which to split the density plot
#' @param yScale vertical scale at which to plot the density. For example, a call with 'yScale = 1' will
#' produce a density curve scaled between 0 and 1  
#' @param fillBackground binary specification of whether to fill in the background the outerCol color
#' @return none. Graph is plotted to the current device
#' @author Jason Waddell
#' @export
#' @examples
#' library(RColorBrewer)
#' data <- rnorm(1000)
#' x <- density(data)$x
#' y <- density(data)$y
#' colVec <- brewer.pal(9, "Blues")[3:8]
#' outerCol <- brewer.pal(9, "Blues")[9]
#'
#' oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
#' drawSplitDensity(x, y, colVec = colVec, split = c(-8), 
#'		outerCol = outerCol,
#'		yScale = 0.95, yshift = 0)
drawSplitDensity <- function(x = NULL, y = NULL, densityObj = NULL, yshift = 0, colVec, outerCol,
		lwd = 2, split = NULL, yScale = NULL, fillBackground = FALSE){
	
	if(!is.null(x) && !is.null(y) && !is.null(densityObj))
		warning("'density', 'x', and 'y' were all provided. Using 'density' to draw plot.")
	
	if(length(split) > (length(colVec) - 1))
		warning(paste("length(split) > length(colVec) - 1. Only the lowest ", length(colVec)-1, 
						" split values will be drawn.", sep = ""))
	
	if(!is.null(densityObj)){
		x <- densityObj$x
		y <- densityObj$y
	}
	
	if(!is.null(split)) split <- sort(split)
	x <- as.numeric(x); y <- as.numeric(y)
	if(!is.null(yScale))
		y <- yScale*y/max(y)
	
	if(fillBackground)
		drawPoly(x, y, yshift = yshift, colShade = outerCol, colOut = outerCol)
	if(is.null(split)){
		drawPoly(x, y, colShade = colVec[1], colOut = outerCol, lwd = lwd, yshift = yshift)
	} else if(!(all(x < split[1]) | all(x > split[length(split)])) ){
		ind <- which(x < split[1])
		drawPolyS1(x[ind], y[ind], colShade = colVec[1], colOut = outerCol, 
				lwd = lwd, yshift = yshift, split = split[1])
		
		for(i in 2:(length(split)+1)) {
			
			if(i == length(split)+1){
				ind2 <- which(x >= split[length(split)])
				drawPolyS2(x[ind2], y[ind2], colShade = colVec[i], colOut = outerCol, 
						lwd = lwd, yshift = yshift, split = split[length(split)])
			} else {
				index <- which(x >= split[i-1] & x < split[i])
				drawPolyMiddle(x[index], y[index], 
						colShade = colVec[i], colOut = outerCol, 
						lwd = lwd, yshift = yshift, split = split, i)
			}
		}
	} else if (all(x <= min(split))){
		warnings("All x values are lower than min(split).")
		drawPoly(x, y, colShade = colVec[1], colOut = outerCol, lwd = lwd, yshift = yshift)
	} else {
		drawPoly(x, y, colShade = colVec[1], colOut = outerCol, lwd = lwd, yshift = yshift)
	}
	lines(x, y+yshift, col = outerCol, lwd = lwd)
}








#


