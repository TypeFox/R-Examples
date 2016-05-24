simer<-function(n, r){
	x <- runif(n, 0, n)
	L <- chol(exp(-abs(outer(x, x, FUN="-"))/r))
	e <- as.vector(crossprod(rnorm(n), L))
	data.frame(x=x, y=e)
}

