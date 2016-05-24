solve.pentadiag <- function(a, b,...){
	a <- as.matrix(a)
	if(nrow(a) != ncol(a))
		stop("'a' is not a square matrix!")
	b <- as.vector(b)
	if(length(b) != ncol(a))
		stop("'a' and 'b' should be of same length!")
	n <- length(b)
	aaa <- diag(a)
	bb <- c(a[cbind(1:(n-1), 2:n)], 0)
	c <- c(a[cbind(1:(n-2), 3:n)], 0, 0)
	d <- c(0, a[cbind(2:n, 1:(n-1))])
	e <- c(0, 0, a[cbind(3:n, 1:(n-2))])
	u <- v <- w <- p <- q <- r <- s <- t <- aa <- dd <- rep_len(0,n+2)
	z <- rep_len(0,n)
	BB <- .C("penta", as.integer(n), as.double(aaa), as.double(bb), as.double(c), as.double(d),
	as.double(e), as.double(b), as.double(u), as.double(v), as.double(w), as.double(p), as.double(q),
	as.double(r), as.double(s), as.double(t), as.double(aa), as.double(dd), as.double(z), 
	PACKAGE = "simest")
	z <- as.vector(BB[[length(BB)]])
	return(z)
}