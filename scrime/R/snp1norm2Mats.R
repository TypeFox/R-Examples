snp1norm2Mats<-function(x1,x2,n.cat=NULL){
	if(is.null(n.cat))
		n.cat<-max(x1,na.rm=TRUE)
	if(n.cat!=3)
		stop("Only available for variables with three levels.")
	if(any(is.na(x1)))
		stop("No missing values allowed in n1.")
	listX2<-getListIdentity(x2,n.cat)
	mat<-matrix(0,nrow(x1),nrow(x2))
	for(i in 1:3){
		idX1<-x1==i
		for(j in (1:3)[-i]){
			tmp<-idX1%*%t(listX2[[j]])
			mat<-mat+abs(i-j)*tmp
		}
	}
	rs<-rowSums(!is.na(x2))
	n<-rep(1,nrow(x1))%*%t(rs)
	mat/(2*n)
}


