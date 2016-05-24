coreCA <-
function(DATA,masses=NULL,weights=NULL,hellinger=FALSE,symmetric=TRUE,decomp.approach='svd',k=0){

	DATA_dimensions = dim(DATA)
	
	
	###PERHAPS ALL OF THIS SHOULD OCCUR IN THE CA FUNCTION?
	mRP<-makeRowProfiles(DATA,weights=weights,masses=masses,hellinger=hellinger)
	#rowCenter <- mRP$rowCenter
	#rowProfiles <- mRP$rowProfiles
	#deviations <- mRP$deviations
	#masses <- mRP$masses
	#weights <- mRP$weights

	#pdq_results <- genPDQ(M=masses,deviations,W=weights,decomp.approach=decomp.approach,k=k)
	pdq_results <- genPDQ(datain=mRP$deviations,M=mRP$masses,W=mRP$weights,is.mds=FALSE,decomp.approach=decomp.approach,k=k)	

	
	#Rows, F
	fi <- pdq_results$p * matrix(pdq_results$Dv,nrow(pdq_results$p),ncol(pdq_results$p),byrow=TRUE)
	rownames(fi) <- rownames(DATA)	
	di <- rowSums(fi^2)
	ri <- repmat((1/di),1,pdq_results$ng) * (fi^2)
	ri <- replace(ri,is.nan(ri),0)	
	ci <- repmat(mRP$masses,1,pdq_results$ng) * (fi^2)/repmat(t(pdq_results$Dv^2),DATA_dimensions[1],1)
	ci <- replace(ci,is.nan(ci),0)
	di <- as.matrix(di)

	###this could be cleaned up. But, after I overhaul CA on the whole.
	fj <- repmat(mRP$weights,1,pdq_results$ng) * pdq_results$q * matrix(pdq_results$Dv,nrow(pdq_results$q),ncol(pdq_results$q),byrow=TRUE)
	rownames(fj) <- colnames(DATA)		
	if(hellinger){
		cj <- (fj^2)/t(repmat(colSums(fj^2),1,nrow(fj)))
	}else{	
		cj <- repmat(mRP$rowCenter,1,pdq_results$ng) * (fj^2)/repmat(t(pdq_results$Dv^2),DATA_dimensions[2],1)
	}
	if(!symmetric){
		fj <- fj * matrix(pdq_results$Dv^-1,nrow(pdq_results$q),ncol(pdq_results$q),byrow=TRUE)
	}
	cj <- replace(cj,is.nan(cj),0)	
	dj <- rowSums(fj^2)
	rj <- repmat((1/dj),1,pdq_results$ng) * (fj^2)
	rj <- replace(rj,is.nan(rj),0)
	dj <- as.matrix(dj)			
	
	return(list(fi=fi,di=di,ci=ci,ri=ri,fj=fj,cj=cj,rj=rj,dj=dj,t=pdq_results$tau,eigs=pdq_results$Dv^2,M=mRP$masses,W=mRP$weights,c= mRP$rowCenter,pdq=pdq_results,X=mRP$deviations,hellinger=hellinger,symmetric=symmetric))
}
