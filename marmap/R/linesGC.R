linesGC <- function(start.points, end.points, n = 10, antimeridian = FALSE, ...) {
	
	if (is.data.frame(start.points)) start.points <- as.matrix(start.points)
	if (is.data.frame(end.points)) end.points <- as.matrix(end.points)
		
	if (is.vector(start.points)) {
		if (length(start.points) != 2) {
			stop("'start.points' must be a vector of length 2 or a 2-column matrix or data.frame")
			} else {
				start.points <- t(as.matrix(start.points))		
			}

	}

	if (is.vector(end.points)) {
		if (length(end.points) != 2) {
			stop("'end.points' must be a vector of length 2 or a 2-column matrix or data.frame")
			} else {
				end.points <- t(as.matrix(end.points))		
			}

	}

	if (!is.matrix(start.points)) stop("'start.points' must be a vector of length 2 or a 2-column matrix or data.frame")
	if (!is.matrix(end.points)) stop("'end.points' must be a vector of length 2 or a 2-column matrix or data.frame")
	
	for (i in 1:nrow(start.points)) {

		# Generate intermediate points along the path
		ll <- suppressWarnings(geosphere::gcIntermediate(start.points[i,], end.points[i,], n=n, addStartEnd=TRUE))
		# print(ll)

		# Take care of antimeridian crossing
		if (antimeridian) {
			ll[ll[,1]<0,1] <- ll[ll[,1]<0,1] + 360
		}
	
		lines(ll, ...)
	}
	
}
