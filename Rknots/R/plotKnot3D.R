# Script comments and history
# 2011
# Feb 11, 2011
# 2:00:51 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

plotKnot3D <- function (points3D, ends = c(), text = FALSE, showNC = FALSE, colors = list(), ...) {
	n <- nrow(points3D)
	ncomp <- length(ends) + 1
	if( identical(colors, list()) ) 
		colors <- as.list( 1 : ncomp) #not supplied by the user
	exends <- c(0, ends, n)
	intervals <- list()
	for (k in 1 : ncomp) 
		intervals[[k]] <- (exends[k] + 1) : (exends[k + 1])
	for (i in 1 : length(intervals)) 
		plotComponent(points3D[intervals[[i]], ], text, showNC,
				col = colors[[i]],  ...)
}

plotComponent <- function (points3D, text = FALSE, showNC = FALSE, ...) {
	if (missing(points3D)) 
		stop("Argument 'points3D' missing, with no default\n")
	if (!is.logical(text) | !is.logical(showNC))
		stop("Argument 'text' and 'showNC' must be logical\n")
	lines3d(points3D, ...)
	material3d(color = c("orange", "brown"), lit = TRUE, 
			point_antialias = TRUE, line_antialias = TRUE)
	bg3d(sphere = TRUE, back = "filled", color= "#887777")
	#bg3d(sphere = FALSE, back = "filled", color="white")
	lines3d(points3D, ...)
	spheres3d(points3D, ...)
	if (text) {
		texts3d(points3D[2 : (nrow(points3D) - 1), ], 
				texts = as.character(2 : (nrow(points3D) - 1)),
				cex = 0.8, adj = 1.35, col = "black")
	}
	if (showNC)
		texts3d(points3D[c(1, nrow(points3D)), ], texts = c("N", "C"),
				cex = 1.5, adj = 1.35, col = "red")
}
