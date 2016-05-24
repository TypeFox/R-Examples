ternaryDiagLines <- function(x, ...){
	s <- rowSums(x)
	if (any(s <= 0)) 
		stop("rowSums of the input data x must be positive.")
	x <- x/s
	top <- sqrt(3)/2
	xp <- x[, 2] + x[, 3]/2
	yp <- x[, 3] * top
	lines(xp, yp, ...)
}
