subjETSDIF <- function(x){

	ETSs <- list(length(x$DIF))
	
	for (i in 1:length(x$DIF)){
		ETSs[[i]] <- subjETS(x$DIF[[i]])
	}
	
	return(ETSs)
	
}