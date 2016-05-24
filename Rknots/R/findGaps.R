##Script generated in:
# 2011
# 9:19:52 AM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

findGaps <- function(points3D, cutoff = 7) {
	
	deuc <- function(x,y) 
		sum((x - y)^2)^(1/2)
	
	n <- nrow(points3D)
	d.vector <- c()
	for (i in 1 : (n - 1))
		d.vector[i] <- deuc(points3D[i, ], points3D[i + 1, ])
	cat("Summary of the distance vector \n")
	print(summary(d.vector))
	gaps <- which(d.vector >= cutoff)
	if(!identical(gaps, integer(0)))
		cat("Chain split at position \n", gaps, "\n")
	else
		cat("No gap found \n")
	split.points3D <- list()
	breaks <- c(0, gaps, n)
	for (i in 1 : (length(breaks) - 1))
		split.points3D[[i]] <- 
				points3D[(breaks[i] + 1) : breaks[i + 1], ]
	return(split.points3D)
}
