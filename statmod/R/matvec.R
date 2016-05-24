matvec <- function(M,v) {
#	Multiply the columns of matrix by the elements of a vector,
#	i.e., compute M %*% diag(v)
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[2]) stop("matvec: Dimensions do not match")
	t(v * t(M))
}

vecmat <- function(v,M) {
#	Multiply the rows of matrix by the elements of a vector,
#	i.e., compute diag(v) %*% M
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[1]) stop("vecmat: Dimensions do not match")
	v * M
}
