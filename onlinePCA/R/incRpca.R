incRpca <- function (lambda, U, x, n, f = 1/n, q = length(lambda), center, tol = 1e-07) 
{
    k <- length(lambda)
	q <- as.integer(q)
    if (ncol(U) != k) 
        stop("length(lambda) != ncol(U)")
    if (nrow(U) != length(x)) 
        stop("length(x) != nrow(U)")
    if (!missing(center)) 
        x <- x - center
    if (missing(tol))
		tol <- sqrt(.Machine$double.eps)
		
   	lambda <- (1-f) * lambda
   	x <- sqrt(f) * x
    xhat <- crossprod(U, x)
    x <- x - U %*% xhat
    normx <- sqrt(sum(x^2))
    if (normx >= tol) {
    	k <- k+1L
    	lambda[k] <- 0
    	xhat[k] <- normx
    	U <- c(U,x/normx)
    	dim(U) <- c(length(x),k)
    }

    eig <- eigen(diag(lambda) + tcrossprod(xhat), TRUE)

    if (q<k) {
    	length(eig$values) <- q	
    	eig$vectors <- eig$vectors[,1:q]
    	}

    return(list(values = eig$values, vectors = U %*% eig$vectors))
}
