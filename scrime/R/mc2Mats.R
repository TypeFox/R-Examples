`mc2Mats` <-
function(x1,x2,use.smc=TRUE,n.cat=NULL){
	if(is.null(n.cat))
		n.cat<-max(x1,na.rm=TRUE)
	listX1<-getListIdentity(x1,n.cat)
	listX2<-getListIdentity(x2,n.cat)
	mat.dist<-matrix(0,nrow(x1),nrow(x2))
	for(i in 1:n.cat){
		tmp<-listX1[[i]]%*%t(listX2[[i]])
		mat.dist<-mat.dist+tmp
	}
	nona1<-!is.na(x1)
	nona2<-!is.na(x2)
	n<-nona1%*%t(nona2)
	if(use.smc)
		return(1-mat.dist/n)
	mat.exp<-matrix(0,nrow(x1),nrow(x2))
	for(i in 1:n.cat){
		tmp1<-listX1[[i]]%*%t(nona2)
		tmp2<-nona1%*%t(listX2[[i]])
		tmp<-tmp1*tmp2/n
		mat.exp<-mat.exp+tmp
	}
	(n-mat.dist)/(n-mat.exp)
}

