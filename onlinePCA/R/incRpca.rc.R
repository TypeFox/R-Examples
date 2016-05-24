incRpca.rc <- function (lambda, Ut, Us, x, n, f = 1/n, center, tol = 1e-7) 
{
	q <- length(lambda)
    if (ncol(Ut) != q) 
        stop("length(lambda) != ncol(Ut)")
    if (nrow(Us) != q || ncol(Us) != q)
    	stop("nrow(Us) != q or ncol(Us) != q")
    if (nrow(Ut) != length(x)) 
        stop("length(x) != nrow(Ut)")
	if (!missing(center))
		x <- x - center 

	x <- sqrt(f) * x
    lambda <- (1-f) * lambda
    xhat <- crossprod(Us, crossprod(Ut,x))
    x <- x - Ut %*% (Us %*% xhat)
    normx <- sqrt(sum(x^2))

    if (normx < tol) {
    	eig <- eigen(diag(lambda) + tcrossprod(xhat), TRUE)
    	return(list(values = eig$values, Ut = Ut, Us = Us %*% eig$vectors))
    }  
        	   	 
    lambda[q+1L] <- 0
    xhat[q+1L] <- normx
    eig <- eigen(diag(lambda) + tcrossprod(xhat), TRUE)
	eig$vectors[1:q,] <- Us %*% eig$vectors[1:q,]
	Ut <- c(Ut,x/normx)
	dim(Ut) <- c(length(x), q+1L)
	return(list(values = eig$values, Ut = Ut, Us = eig$vectors))
}

