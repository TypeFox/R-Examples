#
# given a dataset and a threshold for 'replicate similarity'
#
# returns list of all appropriate groups of (non-identical) replicate samples
# (for use in calculating error model parameters)
###################################################################

findReplicateGroups <- function(distMatrix,theta=0.05){

	passingRows <- which(apply(distMatrix,MARGIN=1,min,na.rm=TRUE)<theta)
	replicateGroups <- list()
	while(length(passingRows)>0){
		a <- passingRows[1]
		theseDists <- distMatrix[a,]
		B <- intersect(which(theseDists>0),intersect(which(!is.na(theseDists)),which(theseDists<qnorm(theta,mean=mean(distMatrix[upper.tri(distMatrix)]),sd=sd(distMatrix[upper.tri(distMatrix)])))))
		if(length(B)>0){
			replicateGroups[[length(replicateGroups)+1]] <- unique(c(a,B))
		}
		passingRows <- setdiff(passingRows,c(a,B))
	}

	# if no groups, include all samples
	if(length(replicateGroups)==0){
		replicateGroups[[1]] <- 1:nrow(distMatrix)
	}

	replicateGroups
}
