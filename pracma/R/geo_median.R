##
##  g e o _ m e d i a n . R  Geometrical Median
##


geo_median <- function(P, tol = 1e-07, maxiter = 200) {
	stopifnot(is.numeric(P))
	if (!is.matrix(P))
	    stop("Argument 'P' must be a matrix (of points in R^n).")

	m <- nrow(P); n <- ncol(P)
    if (n == 1)
        return(list(p = median(P), d = sum(abs(P - median(P))),
                    reltol = 0, niter = 0))

	p0 <- apply(P, 2, mean)
	p1 <- p0 + 1

    iter <- 1
	while(max(abs(p0 - p1)) > tol && iter < maxiter) {
		iter <- iter + 1
		p0 <- p1
	    s1 <- s2 <- 0
	    for (j in 1:m) {
		    d <- Norm(P[j, ] - p0)
		    s1 <- s1 + P[j, ]/d
		    s2 <- s2 + 1/d
	    }
	    p1 <- s1 / s2
	}
	if (iter >= maxiter)
        warning("Maximum number of iterations reached; may not converge.")

    d <- 0
    for (j in 1:m)
        d <- d + Norm(P[j, ] - p1)

	return( list(p = p1, d = d, reltol = max(abs(p0 - p1)), niter = iter) )
}
