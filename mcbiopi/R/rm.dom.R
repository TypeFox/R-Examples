rm.dom<-function(mat,col=FALSE,dom=TRUE){
	if(col)
		mat<-t(mat)
	mat<-mat[!duplicated(mat),,drop=FALSE]
	row.dom<-mat%*%t(mat)==rowSums(mat)
	ids.row<-if(dom) colSums(row.dom)==1 else rowSums(row.dom)==1
	mat<-mat[ids.row,,drop=FALSE]
	#if(sum(ids.row)==1)
	#	mat<-matrix(mat,ncol=1,dimnames=list(rownames(mat),colnames(mat)[ids.row]))
	if(col)
		mat<-t(mat)
	mat
}

