S_conditional <- function(t, z, beta, lambda, a){
     
	result <- exp((-1) * exp(z %*% beta) * lambda * (t ^ a))
     return(result)
}









#