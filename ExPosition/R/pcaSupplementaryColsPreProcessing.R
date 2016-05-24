pcaSupplementaryColsPreProcessing <- function(SUP.DATA=NULL,center=TRUE,scale=TRUE,M=NULL){
	
	if(is.null(SUP.DATA)){
		stop('Must provide supplemental data')
	}
	if(is.null(M)){
		M <- rep(1,nrow(SUP.DATA)) #you need to choose or else I make nothing happen...
	}else if(length(M)!=nrow(SUP.DATA)){
		stop('Length of M does not match row dim of SUP.DATA')
	}
	#return(t(expo.scale(SUP.DATA,center=center,scale=scale)) * t(repmat(M,1,ncol(SUP.DATA))))
	return(t(expo.scale(SUP.DATA,center=center,scale=scale) * repmat(M,1,ncol(SUP.DATA))))
	
}