hybridHclust <- function(x,themc=NULL,trace=FALSE){
	if (trace & is.null(themc)) cat('finding mutual clusters\n')
	if (is.null(themc)) themc <- mutualCluster(x)

	if (trace) cat('temporarily fusing together data belonging to mutual clusters\n')
	xfuse <- x
	for (i in 1:length(themc)){
		mu <- apply(xfuse[themc[[i]],],2,mean,na.rm=TRUE)
		for (j in themc[[i]]) xfuse[j,] <- mu
	}

	if (trace) cat('running tsvq with MCs constrained to remain together\n')
	# K = n - (num points belonging to MCs) + (num  of MCs) 
	#   = num distinct rows of xfuse
	# This is done to prevent tsvq from trying to subdivide identical points
	fdat.tsvq <- tsvq(xfuse,
		K=nrow(xfuse) - length(unlist(themc)) + length(themc),
		row.labs=1:nrow(xfuse),as.hclust=FALSE,trace=trace)
	if (trace) cat('redoing tsvq within each mutual cluster\n')
	result <- tsvq2hclust(redo.fused.tsvq(fdat.tsvq,x,trace=trace))
	return(result)
}
