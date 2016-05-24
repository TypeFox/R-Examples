# Script comments and history
# 2011
# 5:35:25 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

msr <-
function (points3D, ends = c(), n = 100) 
{
	if (missing(points3D)) 
		stop("msr: Argument 'points3D' missing, with no default\n")
	if (!is.numeric(points3D)) 
		stop("msr: Argument 'points3D' must be numeric\n")
	if (!is.null(ends)) {
		ends <- sort(ends)
		if (!is.numeric(ends)) {
			ends <- as.numeric(ends)	
			warning("msr: Argument 'ends' must be numeric, converted to numeric from string\n")
		}
		if (any(ends <= 0)) {
			stop("msr: Argument 'ends' contains negative or 0 entries\n")
		}
	}
	if (!is.numeric(n)) {
		n <- as.numeric(n)
		warning("msr: Argument 'n' must be numeric, converted to numeric from string\n")
	}
	M <- intersectionMatrix(points3D, ends)
	l <- quote(nrow(points3D))
	k <- length(ends) + 1
	for (i in 1 : n) {
		if (eval(l) <= 2 * k) 
			break
		p <- nrow(points3D)
		for (b in 1 : (eval(l) - 1)) {
			temp <- grm(points3D, b, ends, M)
			points3D <- temp$points3D
			ends <- temp$ends
			M <- temp$M
		}
		if (p == eval(l)) 
			break
	}
	return(list(points3D = points3D, ends = ends, M = M))
}


msrFast <-
		function (points3D, ends = c(), n = 100) 
{
	M <- intersectionMatrix(points3D, ends)
	l <- quote(nrow(points3D))
	k <- length(ends) + 1
	for (i in 1 : n) {
		if (eval(l) <= 2 * k) 
			break
		p <- nrow(points3D)
		for (b in 1 : (eval(l) - 1)) {
			temp <- grm(points3D, b, ends, M)
			points3D <- temp$points3D
			ends <- temp$ends
			M <- temp$M
		}
		if (p == eval(l)) 
			break
	}
	return(list(points3D = points3D, ends = ends, M = as.matrix(M)))
}