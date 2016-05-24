ColorsNames<-function(MatrixColors,cols=NULL){
	Names=matrix(0,nrow=dim(MatrixColors)[1],ncol=dim(MatrixColors)[2])
	for (i in 1:dim(MatrixColors)[1]){
		for(j in 1:dim(MatrixColors)[2]){
			temp=MatrixColors[i,j]	
			Color=cols[temp]
			Names[i,j]=Color
		}
	}
	return(Names)
}
