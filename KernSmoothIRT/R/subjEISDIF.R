subjEISDIF <- function(x){

	EISs <- list(length(x$DIF))
	
	for (i in 1:length(x$DIF)){
		EISs[[i]] <- subjEIS(x$DIF[[i]])
	}
	
	return(EISs)
	
}