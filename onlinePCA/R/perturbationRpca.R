perturbationRpca <- function (lambda, U, x, n, f = 1/n, center, sort = TRUE) 
{
	stopifnot(f >= 0 & f <= 1) 
	q <- ncol(U)
	d <- length(x) 
    stopifnot(length(lambda) == q)
    stopifnot(nrow(U) == d)
    if (!missing(center)) 
    	x <- x - center

    lambda <- (1-f) * lambda
    z <- sqrt(f) * crossprod(U,x)
    z2 <- z * z
	num <- tcrossprod(z)  
	den <- matrix(lambda + z2, q, q, byrow = TRUE) - 
		matrix(z2 + lambda^2, q, q)
	V <- num / den
	diag(V) <- 1
	U <- U %*% V
    sigma2 <- .colSums(U * U, d, q)
	lambda <- (lambda + z2) * sigma2
    U <- U * rep.int(1/sqrt(sigma2), rep.int(d,q))
    
    if (sort) {
    	ind <- order(lambda, decreasing = TRUE)
        if (!identical(ind, 1:q)) {
	    	lambda <- lambda[ind]
    		U <- U[,ind]
    	}
	}
	
	list(values = lambda, vectors = U)
}
