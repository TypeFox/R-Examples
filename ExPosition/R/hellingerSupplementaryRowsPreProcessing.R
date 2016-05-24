hellingerSupplementaryRowsPreProcessing <- function(SUP.DATA,center=NULL){
	
	if(is.null(center)){
		#stop("Hellinger supplementary rows require a center (from, e.g., active data)")
		print('No center for Hellinger. Computing center from SUP.DATA')
		hell.preproc <- makeRowProfiles(SUP.DATA,hellinger=TRUE)$deviations
	}
	else{
		hell.preproc <- rowNorms(SUP.DATA,type="hellinger")
		hell.preproc <- hell.preproc - matrix(center,nrow(SUP.DATA),ncol(SUP.DATA),byrow=TRUE)
	}
	
	#test.fi <- supplementalProjection(hell.preproc,f.scores=test.res$ExPosition.Data$fj,Dv=test.res$ExPosition.Data$pdq$Dv,symmetric=FALSE)$f.out

}