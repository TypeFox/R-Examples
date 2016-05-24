affinity <-
function(data, index, column=NULL, R = NULL, weighted=FALSE){
	tmp <- abs(data - matrix(data[index, ], nrow=nrow(data), ncol=ncol(data), byrow=TRUE)) < 1
	if(is.null(column) | is.null(R)){
		if(weighted)warning("Correlation Matrix and/or Column Number not provided, switching to Unweighted Affinity\n")
	}
	if(!weighted){
		affinity <- rowMeans(tmp, na.rm=TRUE)
	}
	if(weighted){
		tmp[which(is.na(tmp), arr.ind=T)] <- FALSE
		affinity <- tmp %*% R[ ,column]
	}
	affinity[index] <- 0
	affinity
}
