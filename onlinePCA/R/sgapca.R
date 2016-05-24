sgapca <- function (lambda, U, x, gamma, q = length(lambda), center, type = c("exact","nn"), sort = TRUE) 
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

	gamma <- rep_len(gamma,k)
	y <- as.numeric(crossprod(U,x))
	type <- match.arg(type)
	U <- switch(type, exact = sgapca_exC(U,x,y,gamma), 
		nn = sgapca_nnC(U,x,y,gamma))

    if (!missing(lambda)) {
        lambda <- (1 - gamma) * lambda + gamma * y * y
        if (sort) {
            ix <- order(lambda, decreasing = TRUE)
            if (!identical(ix, 1:k)) {
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
