# added 2010-06-22; by M. Friendly, modified by  J. Fox

wcrossprod <- function(x, y, w) {
	if (is.vector(x)) x <- as.matrix(x)
	if (!missing(y)){
		if (is.vector(y)) y <- as.matrix(y)
		if (nrow(x) != nrow(y)) stop("x and y not conformable")
	}
	if (missing(w)) {
		if (missing(y)) return(crossprod(x)) else return(crossprod(x, y))
	}
	else if (length(w)==1 || (is.vector(w) && sd(w) < sqrt(.Machine$double.eps))) {
		if (missing (y)) return(w[1]*crossprod(x)) else return(w[1]*crossprod(x, y))
	}
	else {
		if (is.vector(w)) {
			if (length(w) != nrow(x)) stop("w is the wrong length")
			if (missing(y)) return(crossprod(x, w*x)) else return(crossprod(x,  w*y))
		}
		else {
			if (nrow(w) != ncol(w) || nrow(w) != nrow(x)) stop("w is the wrong dimension")
			if (missing(y)) return(crossprod(x, w %*% x)) else return(crossprod(x, w %*% y))
		}
	}
}
