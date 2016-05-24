caNorm <-
function(X,X_dimensions,colTotal,rowTotal,grandTotal,weights=NULL,masses=NULL){

	rowCenter = colTotal/grandTotal
	if(is.null(masses)){
		masses = rowTotal/grandTotal
	}
	#M = diag(masses)
	if(is.null(weights)){
		weights = rowCenter^-1
	}
	#W = diag(weights) 
	#rowProfiles <- X/((rowSums(X))%*%matrix(1,1,ncol(X)))
	rowProfiles <- rowNorms(X,type='ca')
	#deviations = rowProfiles - (repmat(1,X_dimensions[1],1)%*%rowCenter)
	deviations <- rowProfiles - matrix(rowCenter,X_dimensions[1],X_dimensions[2],byrow=TRUE)
	#return(list(rowCenter=rowCenter,masses=masses,M=M,weights=weights,W=W,rowProfiles=rowProfiles,deviations=deviations))
	return(list(rowCenter=rowCenter,masses=masses,weights=weights,rowProfiles=rowProfiles,deviations=deviations))
}
