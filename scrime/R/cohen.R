`cohen` <-
function(x,dist=FALSE){
	out<-computeContCells(x,justDiag=TRUE)
	aR<-rowSums(out$mat.obs)
	eR<-rowSums(out$mat.exp)
	if(any(is.na(x))){
		naIdentifier<-!is.na(x)
		n<-naIdentifier%*%t(naIdentifier)
		n<-n[lower.tri(n)]
	}
	else
		n<-ncol(x)
	coef<-(aR-eR)/(n-eR)
	mat<-matrix(0,nrow(x),nrow(x))
	mat[lower.tri(mat)]<-coef
	mat<-mat+t(mat)
	diag(mat)<-1	
	rownames(mat)<-colnames(mat)<-rownames(x)
	if(dist)
		mat<-1-mat
	mat
}

