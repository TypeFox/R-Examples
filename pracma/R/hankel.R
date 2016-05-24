##
##  h a n k e l . R
##


hankel <- function(a, b) {
    if (!is.vector(a))
        stop("Argument 'a' must be a numeric or complex vector.")
    n <- length(a)
	if (missing(b)) b <- c(a[n], rep(0, n-1))
    if (!is.vector(b))
        stop("Argument 'b' must be a numeric or complex vector.")
    m <- length(b)
	
	if (a[n] != b[1])
	    warning("a[n] not equal to b[1], b[1] set to a[n].")

	H <- matrix(0, n, m)
	for (i in 2:(n+m))
		H[row(H)+col(H) == i] <- if (i <= n+1) a[i-1] else b[i-n]
	return(H)
}
