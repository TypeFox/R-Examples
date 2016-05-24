ccipca <- function(lambda, U, x, n, q=length(lambda), 
	l=2, center, tol = 1e-8, sort = TRUE)
{
	n <- as.integer(n)
	q <- min(q, n+1L, length(lambda)+1L)
    stopifnot(ncol(U) == length(lambda)) 
    stopifnot(nrow(U) == length(x)) 
    stopifnot(l >= 0 && l <= n)
    if (!missing(center)) 
        x <- x - center
	result <- ccipca_C(lambda, U, x, n, q, l, tol)

	if (sort) {
		ix <- order(result$values, decreasing=TRUE)
		if (!identical(ix, 1:q)) {
			result$values <- result$values[ix]
		result$vectors <- result$vectors[,ix]
		}
	}
	
	return(result)
}

