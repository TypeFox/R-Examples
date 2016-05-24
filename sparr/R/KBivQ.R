
KBivQ <- function(X,type="spher"){

	if(!is.data.frame(X) && !is.vector(X)) stop("X must be a numeric vector of length 2 or a 2-column data frame")
	
	if(is.vector(X)&&length(X)>2){
		X <- X[1:2]
		warning("X was a vector with > 2 components... truncated to the first two")
	}
	
	if(is.data.frame(X)&&ncol(X)>2){
		X <- X[,1:2]
		warning("X was a data frame with > 2 columns... truncated to the first two")
	}

    if(is.data.frame(X)){
        if(type=="spher"){
            u <- sqrt(X[,1]^2+X[,2]^2)
            result <- (3/pi)*(16/15)*(15/16)*(1-u^2)^2 * (abs(u)<=1)
        } else {
            u1 <- X[,1]
            u2 <- X[,2]  
            result <- ((15/16)*(1-u1^2)^2)*((15/16)*(1-u2^2)^2) * (abs(u1)<=1) * (abs(u2)<=1)
        }
    } else {
        if(type=="spher"){
            u <- sqrt(X[1]^2+X[2]^2)
            result <- (3/pi)*(16/15)*(15/16)*(1-u^2)^2 * (abs(u)<=1)
        } else {
            u1 <- X[1]
            u2 <- X[2]  
            result <- ((15/16)*(1-u1^2)^2)*((15/16)*(1-u2^2)^2) * (abs(u1)<=1) * (abs(u2)<=1)
        }
    }
    return(result)
}
