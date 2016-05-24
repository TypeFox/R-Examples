KBivN <- function(X){
	if(!is.data.frame(X) && !is.vector(X)) stop("X must be a numeric vector of length 2 or a 2-column data frame")
	
	if(is.vector(X)&&length(X)>2){
		X <- X[1:2]
		warning("X was a vector with > 2 components... truncated to the first two")
	}
	
	if(is.data.frame(X)&&ncol(X)>2){
		X <- X[,1:2]
		warning("X was a data frame with > 2 columns... truncated to the first two")
	}
	
    if(is.data.frame(X)) result <- (1/(2*pi))*exp(-.5*(X[,1]^2+X[,2]^2))
    else result <- (1/(2*pi))*exp(-.5*(X[1]^2+X[2]^2))
    return(result)
}
