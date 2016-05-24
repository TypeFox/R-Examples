# Returns a list of connected data sets as numeric data frames D 
# with three columns (x, factor 1, factor 2) or a n*m matrix M,
# where the n rows correspond to n levels of factor 2 and m columns
# correspond to m levels of factor the respective factors.
# Output as data frame or matrix, depending on input
getConnectedSets <- function(M) {
	if(!is.data.frame(M) && !is.matrix(M)) {
		stop("getConnectedSets: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(M) && length(M)!=3) {
		stop("getConnectedSets: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(M)) {
		f1 <- factor(M[[2]])
		f2 <- factor(M[[3]])
		M <- raw2matrix(M)
		out <- 1
	}
	else { out <- 0 }

	sets <- connectedSets(M) # find connected sets

	csets <- c()
	for(i in sort(unique(sets$colclasses[which(sets$colclasses!=-1)]))) {
		cs <- M[which(sets$rowclasses==i),which(sets$colclasses==i)]
		# convert to row or column matrix, if a vector is returned
		# if only one level of factor f2, return column matrix
		if(length(which(sets$rowclasses==i)) == 1 & length(which(sets$colclasses==i)) > 1) cs <- t(cbind(cs))
		# if only one level of factor f1, return row matrix
		if(length(which(sets$rowclasses==i)) >= 1 & length(which(sets$colclasses==i)) ==1) cs <- cbind(cs)

		if(out == 1) { # return matrix
			cs <- matrix2raw(cs,as.numeric(levels(f1)[which(sets$rowclasses==i)]),as.numeric(levels(f2)[which(sets$colclasses==i)]))
		}
		csets <- c(csets,list(n=cs))
		names(csets)[[i]] <- paste("cs_",i,sep="")
	}
	return(csets)
}
