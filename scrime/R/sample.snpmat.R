`sample.snpmat` <-
function(n.row,n.col,maf,raf,exclude,vec.equal){
	mat.snp<-matrix(NA,2*n.row,n.col)
	for(i in 1:n.col){
		tmp<-sample(0:1,4*n.row,replace=TRUE,prob=c(raf[i],maf[i]))
		mat.snp[,i]<-tmp[1:(2*n.row)]+tmp[(2*n.row+1):(4*n.row)]
	}
	mat.snp<-checkMatSNP(mat.snp,exclude,vec.equal)
	n.row2<-nrow(mat.snp)
	if(n.row2>=n.row)
		ids.obs<-sample(n.row2,n.row)
	else
		ids.obs<-sample(n.row2,n.row,replace=TRUE)
	mat.snp<-mat.snp[ids.obs,]
	mat.snp
}

