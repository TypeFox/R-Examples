hellingerNorm <-
function(X,X_dimensions,colTotal,rowTotal,grandTotal,weights=NULL,masses=NULL){
	
	if(is.null(masses)){
		masses = rowTotal/grandTotal
	}
	#M = diag(masses)
	if(is.null(weights)){
		weights <- c(matrix(1/ncol(X),1,ncol(X)))
	}
	#W = diag(weights)
	#rowProfiles <- (X/repmat(rowTotal,1,ncol(X)))^(1/2)
	rowProfiles <- rowNorms(X,type='hellinger')	
	rowCenter <- c(t(as.matrix(masses)) %*% rowProfiles)
	#deviations <- rowProfiles - (repmat(1,X_dimensions[1],1)%*%rowCenter)
	deviations <- rowProfiles - matrix(rowCenter,X_dimensions[1],X_dimensions[2],byrow=TRUE)
	#return(list(rowCenter=rowCenter,masses=masses,M=M,weights=weights,W=W,rowProfiles=rowProfiles,deviations=deviations))
	return(list(rowCenter=rowCenter,masses=masses,weights=weights,rowProfiles=rowProfiles,deviations=deviations))	
}
