"make.distmat" <-
function(X) {
	X <- t(X)
	return(sqrt(abs(sweep(sweep(-2*t(X)%*%X, 1, colSums(X*X), "+"), 2, colSums(X*X), "+"))))
}
