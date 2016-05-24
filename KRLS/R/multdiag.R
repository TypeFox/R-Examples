#Function to multiply a square matrix, X, with a diagonal matrix, diag(d)
multdiag <- function(X,d){	
	R=matrix(NA,nrow=dim(X)[1],ncol=dim(X)[2])		
	for (i in 1:dim(X)[2]){
		R[,i]=X[,i]*d[i]	
	}
	return(R)
} 