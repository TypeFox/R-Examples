`getBinGenotype` <-
function(mat){
	n.row<-nrow(mat)
	binmat<-matrix(0,3*n.row,ncol(mat))
	binmat[3*(1:n.row)-2,]<-1*(mat==1)
	binmat[3*(1:n.row)-1,]<-1*(mat==2)
	binmat[3*(1:n.row),]<-1*(mat==3)
	rownames(binmat)<-paste(rep(rownames(mat),e=3),rep(1:3,n.row),sep="_")
	colnames(binmat)<-colnames(mat)
	binmat
}

