`getBinDomRec` <-
function(mat){
	n.row<-nrow(mat)
	binmat<-matrix(0,2*n.row,ncol(mat))
	binmat[2*(1:n.row)-1,]<-1*(mat!=1)
	binmat[2*(1:n.row),]<-1*(mat==3)
	rn<-rownames(mat)
	rownames(binmat)<-paste(rep(rn,e=2),rep(1:2,n.row),sep="_")
	colnames(binmat)<-colnames(mat)
	binmat
}

