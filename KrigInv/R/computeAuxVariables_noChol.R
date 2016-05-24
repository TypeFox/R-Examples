computeAuxVariables_noChol <- function(model){
	
	T <- model@T
	
	z <- backsolve(t(T), model@y-model@F%*%as.matrix(model@trend.coef), upper.tri=FALSE) 	
	model@z <- as.numeric(z)
	
	return(model)
}
