##I almost certainly don't need to send the whole set of data along, do I? This is probably a waste. FIX THIS IN FUTURE RELEASE.
boot.samples <- function(DATA,DESIGN=NULL,constrained=FALSE){
	if(constrained && !is.null(DESIGN)){
		boot.index <- vector()
		for(i in 1:ncol(DESIGN)){
			boot.index <- c(boot.index,sample(which(DESIGN[,i]==1),replace=TRUE))
		}
	}else{
		boot.index <- sample(nrow(DATA),nrow(DATA),TRUE)
	}
	
	return(boot.index)
}