enforceDomainOnParameters <- function(paramsArray, domain) {
# A function to check the if the generated population is within the domain
# if it is not, the parameter that is not respecting the domain will be placed
# inside the closest boundary.
	for(i in 1:length(paramsArray[,1])) {
		for(j in 1:length(paramsArray[1,])) {
			if(paramsArray[i,j] < domain[j,1]) { paramsArray[i,j] <- domain[j,1] }
			if(paramsArray[i,j] > domain[j,2]) { paramsArray[i,j] <- domain[j,2] }
		}
	}
	return(paramsArray)
}