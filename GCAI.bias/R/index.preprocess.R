index.preprocess <-
function(mat.seq,word){

	mat.index<-matrix(0,nrow=nrow(mat.seq),ncol=word*4)

	for(t in 1:nrow(mat.seq)){
		mat.index[t,unlist(mat.seq[t,4+(1:word)])]<-1
	}

	return(mat.index)
}
