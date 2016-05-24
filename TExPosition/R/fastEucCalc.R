fastEucCalc <-
function(x,c){
	
	if(ncol(x) == 1){
		return((x^2) %*% matrix(1,1,nrow(c)) + matrix(1,nrow(x),1) %*% t((c^2))-2 * x %*% t(c))
	}
	else{ 
		x2=colSums(t(x)^2)
		c2=colSums(t(c)^2)
		return(x2 %*% matrix(1,1,nrow(c))+ matrix(1,nrow(x),1) %*% c2-(2 * x %*% t(c)))
	}
}
