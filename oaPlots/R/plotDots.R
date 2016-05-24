#' Adds Points on a Pre-existing Plot using Shifted Locations
#' @param vec numeric vector
#' @param xLeft left x boundary of the point plotting region
#' @param xRight right x boundary of the point plotting region
#' @param ... further arguments to be handed to the points function
#' @return points are added to the current graphics device
#' @author Jason Waddell
#' @examples x <- sample(1:5, size = 25, replace = TRUE)
#' plot(x = -1, y = -1, xlim = c(0.5,1.5), ylim = range(x), 
#'     ylab = "", xlab = "", xaxt = "n")
#' colVec <- c(rep("olivedrab", 15), rep("goldenrod", 5), rep("red", 5))
#' plotDots(vec = x, xLeft = 0.8, xRight = 1.2, pch = 19, 
#'     col = colVec, cex = 2)
#' @export
plotDots <- function(vec = NULL, xLeft = 0.8, xRight = 1.2, ...){
	center <- (xLeft+xRight)/2
	span <- abs(xRight - xLeft)
	
	maxRep <- max(table(vec))
	space <- span/(maxRep-1)
	
	xLocations <- numeric(length = length(vec))
	
	for(i in 1:length(table(vec))){
		numRep <- table(vec)[i]
		tempLoc <- findLocations(n = numRep, space = space, center = center)
		
		tempIndex <- which(vec == sort(unique(vec))[i])
		
		for(j in seq_along(tempIndex))
			xLocations[tempIndex[j]] <- sort(tempLoc)[j]	
			
	}
	
	points(x = xLocations, y = vec,  ...)
}
