`pcc` <-
function(x,dist=FALSE,corrected=TRUE,version=1){
	if(!version%in%(1:3))
		stop("version must be between 1 and 3.")
	if(dist & !corrected)
		stop("corrected must be TRUE, if dist = TRUE.")
	chi2<-rowChisqStats(x, compPval = FALSE)
	chi2<-chi2+t(chi2)
	if(any(is.na(x))){
		naIdentifier<-!is.na(x)
		mat.n<-naIdentifier%*%t(naIdentifier)
	}
	else
		mat.n<-ncol(x)
	pcc<-chi2/(chi2+mat.n)
	mat.rc<-minrc(x)
	mat.rc<-mat.rc/(mat.rc-1)
	diag(pcc)<-if(is.matrix(mat.rc)) 1/diag(mat.rc) else 1/mat.rc
	if(corrected)
		pcc<-mat.rc*pcc
	colnames(pcc)<-rownames(pcc)<-rownames(x)
	if(!dist){
		pcc<-sqrt(pcc)
		return(pcc)
	}
	pcc<-switch(version,sqrt(1-pcc),1-sqrt(pcc),1-pcc)
	pcc
}

