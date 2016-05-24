find_nn_k <-
function(X,k,iLLE=FALSE){
	
	#calculate distance-matrix
	nns <- as.matrix(dist(X))
	
	#get ranks of all entries
	nns <- t(apply(nns,1,rank))
	
	#choose the k+1 largest entries without the first (the data point itself)
	nns <- ( nns<=k+1 & nns>1 )
	
	#optional: improved LLE
	if( iLLE ){
		N <- dim(X)[1]
		n <- dim(X)[2]
		nns2 <- nns
		nns <- as.matrix(dist(X))
		for( i in 1:N){
			if( i%%100==0 ) cat(i,"von",N,"\n")
			for( j in 1:N){
				Mi <- sqrt( sum( rowSums( matrix( c(t(X[i,])) - c(t(X[nns2[i,],])), nrow=k, ncol=n, byrow=TRUE)^2 )^2 ) )
				Mj <- sqrt( sum( rowSums( matrix( c(t(X[j,])) - c(t(X[nns2[j,],])), nrow=k, ncol=n, byrow=TRUE)^2 )^2 ) )
				nns[i,j] <- nns[i,j]/(sqrt(Mi*Mj)*k) 
			}
		}
		nns <- t(apply(nns,1,rank))
		nns <- ( nns<=k+1 & nns>1 )
	}
	
	return (nns)
}

