cv.score <- function(bandwidth, x, y, degree)
{
	spacing <- diff(x)
	if (any(spacing < 0))
		stop("'x' must be increasing") 
	if (!isTRUE(all.equal(min(spacing),max(spacing)))) 
		stop("'x' must be a uniform grid")		
	if (nrow(y) < 2)
		stop("'y' must have at least two rows")
	if (length(x) != ncol(y)) 
		stop("length(x) and ncol(y) must be equal")

	n <- nrow(y)
	N <- ncol(y)
	y.hat  <- apply(y, 1, function(z) locpoly (x = x, y = z, 
			  		bandwidth = bandwidth, gridsize = N, degree = degree)$y)
	mu.hat <- rowMeans(y.hat)
	residuals <- (n/(n-1)) * mu.hat - y.hat / (n-1) - t(y)
	mean(residuals^2)
}
