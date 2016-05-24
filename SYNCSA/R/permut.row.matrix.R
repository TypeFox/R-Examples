permut.row.matrix <-
function(matrix){
	row<-dim(matrix)[1]
	col<-dim(matrix)[2]
	samp<-sample(1:row,row)
	permut.matrix <- matrix[samp,1:col]
	return(permut.matrix)
	}

