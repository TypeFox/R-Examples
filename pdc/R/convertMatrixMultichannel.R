

convert.matrix.multichannel <-
function(X, m, t) {
	warning("Calls to convert.matrix.multichannel(...) is deprecated!")
	return(convertMatrixMultichannel(X, m, t));
}

convertMatrixMultichannel <-
function(X, m, t) {
result <- array(0, dim=c( dim(X)[2]*dim(X)[3] , factorial(m)))

c<-1
for (i in 1:dim(X)[2]) {
	for ( j in 1:dim(X)[3]) {
		result[c, ] <-codebook(X[,i,j],m,t)
		c<-c+1
	}
}

return(result)
}
