rvCoeff <- function(S,T,type=-1){
##private function
matrixTrace <- function(squareMatrix){
	return(sum(diag(squareMatrix)))
}	
	
	rv = 0
	if(sum(dim(S) == dim(T)) != 2){
		print('Dimensions do not match')
		return(NaN)
	}
		

	if(type == 0){
		rv=matrixTrace(t(S) %*% T) / sqrt( matrixTrace(t(S) %*% S) %*% matrixTrace(t(T) %*% T) ) 
	}else if(type == 1){
		rv = t(as.vector(S)) %*% as.vector(T)/sqrt( ((as.vector(S)) %*% as.vector(S)) %*% ((as.vector(T)) %*% as.vector(T)) )
	}else{
		#this is the worst way to do this.
		top = 0
		bottom_left = 0
		bottom_right = 0
		for(i in 1:dim(S)[1]){
			for(j in 1:dim(T)[2]){
				top = top + (S[i,j] * T[i,j])
				bottom_left = bottom_left + S[i,j]^2
				bottom_right = bottom_right + T[i,j]^2
			}
		}
		rv = top/sqrt(bottom_left * bottom_right)
	}
	
	return(as.numeric(rv))

}
