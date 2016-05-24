# Finds connected data sets of a numeric matrix M
# missing values are assumeds to be NA or 0.
# Returns two vectors: 
# rowclasses[0..maxnr-1] : Class number of the respective rows
# colclasses[0..maxnc-1] : Class number of the respective cols
connectedSets <- function(M) {
	if(!is.matrix(M)) stop("connectedSets: first argument must be a matrix. Exiting ...")
	maxnr <- dim(M)[1]
	maxnc <- dim(M)[2]

	# set NA values to 0
	for(i in 1:maxnr) {
		for(j in 1:maxnc) { 
			if(is.na(M[i,j])) { M[i,j] <- 0 }
		}
	}
	res <- .C("Cconnectivity",M=as.vector(t(M),"numeric"),nrows=as.integer(maxnr),ncols=as.integer(maxnc),rowclasses=vector("integer",maxnr),colclasses=vector("integer",maxnc),PACKAGE="pheno")
	
	return(list(rowclasses=res$rowclasses, colclasses=res$colclasses))

}
