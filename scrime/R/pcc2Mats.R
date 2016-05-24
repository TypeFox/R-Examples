`pcc2Mats` <-
function(x1,x2,chisq=FALSE,version=1,n.cat=NULL){
	if(!version%in%(1:3))
		stop("version must be 1, 2, or 3.")
	if(is.null(n.cat))
		n.cat<-max(x1,na.rm=TRUE)
	listX1<-getListIdentity(x1,n.cat)
	listX2<-getListIdentity(x2,n.cat)
	mat<-matrix(0,nrow(x1),nrow(x2))
	nona1<-!is.na(x1)
	nona2<-!is.na(x2)
	n<-nona1%*%t(nona2)
	listX2no1<-vector("list",n.cat)
	for(i in 1:n.cat)
		listX2no1[[i]]<-nona1%*%t(listX2[[i]])
	for(i in 1:n.cat){
		x1<-listX1[[i]]
		x1no2<-x1%*%t(nona2)
		for(j in 1:n.cat){
			Nobs<-x1%*%t(listX2[[j]])
			Nexp<-x1no2*listX2no1[[j]]/n
			Nexp[Nexp==0]<-1
			tmp<-Nobs*Nobs/Nexp
			mat<-mat+tmp
		}
	}
	if(chisq)
		return(mat-n)
	catX1<-compVarLevs(listX1)
	catX2<-compVarLevs(listX2)
	mat.rc<-minrc2Mats(catX1,catX2,n.cat)
	mat.rc<-mat.rc/(mat.rc-1)
	mat<-mat.rc*(mat-n)/mat
	if(any(mat>1))
		mat[mat>1]<-1
	mat<-switch(version,sqrt(1-mat),1-sqrt(mat),1-mat)
	mat
}

