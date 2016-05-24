errorMetric <- function(obs, forec, type="sAPE", statistic="M"){
	
	if( !any(type==c("AE","SE","APE","sAPE")) ) stop("Error in errorMetric function: this error type has not been implemented.")
	if( !any(statistic==c("M","Md","N")) ) stop("Error in errorMetric function: this statistic has not been implemented.")
	if( is.ts(obs) || is.ts(forec) ){
		obs = as.numeric(obs)	
		forec = as.numeric(forec)
	}
	obs = as.matrix(obs)
	forec = as.matrix(forec)
	if(any(dim(obs)!=dim(forec)))  stop("Error in errorMetric function: the dimensions of the vectors are different.")
	
	if(type == "AE")
		errors = abs(obs - forec) 	
	
	if(type == "SE")
		errors = (obs - forec)^2 

	if(type == "APE")
		errors = abs( 100*(obs - forec)/obs )
	
	if(type == "sAPE")
		errors =  abs( 200*(obs - forec)/(abs(obs) + abs(forec)) )	
	
	if(statistic == "M")
		return(  mean(errors, na.rm=TRUE) )
	
	if(statistic == "Md")
		return(  median(errors, na.rm=TRUE) )
	
	return( errors )	

}

