`smc` <-
function(x,dist=FALSE){
	out<-computeContCells(x,computeExp=FALSE,justDiag=TRUE)$mat.obs
	mat<-matrix(0,nrow(x),nrow(x))
	mat[lower.tri(mat)]<-rowSums(out)
	mat<-mat+t(mat)
	if(any(is.na(x))){
		naIdentifier<-!is.na(x)
		n<-naIdentifier%*%t(naIdentifier)
	}
	else
		n<-ncol(x)
	mat<-mat/n
	diag(mat)<-1
	rownames(mat)<-colnames(mat)<-rownames(x)
	if(dist)
		mat<-1-mat
	mat
}

