laplacian <-
function(A, normalised = FALSE){
	
	n = dim(A)[1]
	temp = apply(abs(A), 2, sum)
	D = diag(temp, nrow = n)
	
	temp1 = Reciprocal(sqrt(temp))
	half.D = diag(temp1, nrow = n)
	if(normalised == TRUE) 	return(half.D %*% (D - A) %*% half.D)
	if(normalised == FALSE) return(D - A)
	
}
