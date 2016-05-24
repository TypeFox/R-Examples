test.FastgMCP <- function() {
	# Testing the C interface
	m <- matrix(1/3,nr=4,nc=4)
	diag(m) <- 0
	w <- rep(1/4,4)
	p <- c(0.01, 0.012, 0.07, 0.03)
	a <- 0.05
	result <- gMCP:::fastgMCP(m, w, p, a)
	checkEquals(result$w, c(0,0,0.5,0.5))
	checkEquals(result$m, structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0), .Dim = c(4L, 4L)))
	checkEquals(result$rejected, c(TRUE, TRUE, FALSE, FALSE))
}