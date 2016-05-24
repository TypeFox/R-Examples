pcaSupplementaryRowsPreProcessing <- function(SUP.DATA=NULL,center=TRUE,scale=TRUE,W=NULL){
	if(is.null(SUP.DATA)){
		stop('Must provide supplemental data')
	}
	if(is.null(W)){
		W <- rep(1,ncol(SUP.DATA)) #you need to choose or else I make nothing happen...
	}else if(length(W)!=ncol(SUP.DATA)){
		stop('Length of W does not match column dim of SUP.DATA')
	}
	
	return(expo.scale(SUP.DATA,center=center,scale=scale) * t(repmat(W,1,nrow(SUP.DATA))))
	
}