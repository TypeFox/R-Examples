getEdgeList <- function(m, rm.zero=TRUE) {
	# TODO check if this returns the upper or lower triangle
	
	from_idx <- vector(mode='integer', length=((nrow(m)^2)-nrow(m))/2)
	to_idx <- vector(mode='integer', length=((nrow(m)^2)-nrow(m))/2)
	weight <- vector(mode='numeric', length=((nrow(m)^2)-nrow(m))/2)
	
	result <- .Fortran('getedgelist', mat=as.single(m), nGenes=as.integer(nrow(m)), from_idx=as.integer(from_idx), to_idx=as.integer(to_idx), weight=as.single(weight))
	
	edgeList <- data.frame(From=rownames(m)[result$from_idx], To=colnames(m)[result$to_idx], Weight=as.vector(result$weight))
	
	return(edgeList[edgeList$Weight != 0,])
}
