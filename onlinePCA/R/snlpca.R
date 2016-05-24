snlpca <- function (lambda, U, x, gamma, q = length(lambda), center, type = c("exact", "nn"), sort = TRUE) 
{
	d <- NROW(U)
	k <- NCOL(U)

	stopifnot(length(x)==d)
	if (!missing(lambda))
		stopifnot(length(lambda)==k)		
    if (!missing(center)) 
    	x <- x - center
	if (!is.matrix(U))
		U <- as.matrix(U)

	y <- as.numeric(crossprod(U, x))
	gamma <- rep_len(gamma, k)
    gamy <- gamma * y
	type <- match.arg(type)		

	if (type == "exact") {
	    U <- U + tcrossprod(x, gamy)
	    eig <- eigen(crossprod(U), TRUE)
	    nonzero <- which(eig$values > sqrt(.Machine$double.eps))
	    iS <- eig$vectors[, nonzero] %*% 
	    	(t(eig$vectors[, nonzero])/sqrt(eig$values[nonzero]))
	    U <- U %*% iS
		if (length(nonzero) < ncol(U)) 
        	warning(paste("Matrix 'U' is not full rank. Returning", 
            	length(nonzero), "PC."))
	} else if (type == "nn") {
		U <- U + tcrossprod(x - U %*% y, gamy)
	}

    if (!missing(lambda)) {
        lambda <- (1 - gamma) * lambda + gamma * y * y
        if (sort) {
            ix <- order(lambda, decreasing = TRUE)
            if (!identical(ix, 1:q)) {
                lambda <- lambda[ix]
                U <- U[, ix]
            }
        }
		if (q<k) 
			length(lambda) <- q
    } else lambda <- NULL
 
     if (q<k)
    	U <- U[,1:q]

    return(list(values = lambda, vectors = U))
}
