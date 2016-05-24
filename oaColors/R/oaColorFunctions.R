
#' Generate a Color Wheel
#' @param colVec vector of colors
#' @author Jason Waddell
#' @importFrom grid grid.newpage pushViewport viewport grid.points gpar
#' @importFrom MASS mvrnorm
#' @export
colorWheel <- function(colVec = NULL){
	
	w = 2
	k = length(colVec)
	grid.newpage()
	plotvp <- viewport(x=0.1,
			y=0.1, xscale = c(-w, w), yscale = c(-w, w), 
			width = 0.85, height = 0.85, 
			name="plotRegion", 
			just = c("left", "bottom"))
	pushViewport(plotvp)
	
	sigMat <- matrix(c(0.15, 0, 0, 0.15), ncol = 2)
	
	n = 30
	for(i in 1:k){
		x = cos(2*pi * (i-1)/k )
		y = sin(2*pi * (i-1)/k) 
		
		temp <- mvrnorm(n, mu = c(x,y), Sigma = sigMat)
		grid.points(x = temp[,1], y = temp[,2], pch = 19,
				gp = gpar(col = colVec[i]) )
	}	
}





#' Pick One or More OA Colors
#' @param color a character vector of color names; possible values are "red", "orange", "yellow", "green", "cyan", "blue", "pink", 
#'       "limegreen", "purple", "black", "white", "grey" or "gray" 
#' @param alpha transparency level for the color(s)
#' @return character vector of colors 
#' @author Tobias Verbeke
#' @importFrom grDevices hcl
#' @export
oaColors <- function(color = NULL, alpha = 1.0){

	colPalette1 <- hcl(seq(0, 360, by = 10 ), l = 50, c = 250, alpha = alpha)
	colPalette2 <- hcl(seq(250, 270, by = 1), l = 50, c = 360, alpha = alpha)
	
	colorValues <- c(colPalette1[2], colPalette2[6], colPalette2[15], colPalette1[12], 
			colPalette1[21], colPalette1[27], colPalette1[29], 
			rainbow(60, alpha = alpha)[17], rainbow(360, alpha = alpha)[278], 
			hcl(180, c = 0, l = 0, alpha = alpha), hcl(269, l = 100, c = 90, alpha = alpha), 
			hcl(180, l = 60, c = 0, alpha = alpha), hcl(180, l = 60, c = 0, alpha = alpha)
			)
	names(colorValues) <- c("red", "orange", "yellow", "green", "cyan", "blue", "pink", 
			"limegreen", "purple", "black", "white", "grey", "gray")
	
	out <- colorValues[color]
	
	if(any(is.na(out)))
		stop(
				cat("Error: ", color[which(is.na(out))], " is not in the list of oaColors. Current colors include: ", "\n         ", 
						paste(paste(names(colorValues), collapse = ", "), sep = ""), "\n", sep = "")
		)
		

	return(out)	
}


#' Generate a Palette of OA Colors
#' @param numColors number of colors to be contained in the palette
#' @param alpha transparency level of the colors
#' @return vector of colors
#' @author Jason Waddell
#' @importFrom grDevices rainbow
#' @importFrom RColorBrewer brewer.pal
#' @export
oaPalette <- function(numColors = NULL, alpha = 1.0){
	
	fullPaletteNames <- c("red", "blue", "green", "orange", "pink", "cyan", "yellow", "limegreen", "purple")
	
	if(numColors > 12){
		samplePalette <- rainbow(numColors+1)[1:numColors]
	} else if(numColors > 9){
		samplePalette <- brewer.pal(numColors,"Paired")
	} else {
		samplePaletteNames <- fullPaletteNames[1:numColors]
		samplePalette <- oaColors(samplePaletteNames, alpha = alpha)
	}
	return(samplePalette)
	
}


