dantzig.selector <- function(lambdalist, BETA0, lambda){

	# if BETA0 dimensions do not match, throw error
	if (length(lambdalist) != dim(BETA0)[2]) {
		stop("BETA0 and lambdalist dimensions are incompatible \n")
	}

	if (lambdalist[length(lambdalist)] > lambda) {
		beta0<-BETA0[,length(lambdalist)]
	} 
	else {
		beta0<-BETA0[,which.max(lambdalist <= lambda)]
	}
	return(beta0)	
}