counts.preprocess <-
function(mat.counts){

	counts<-mat.counts[,5:6]
	pos<-mat.counts[,1:4]
	
	return(list(counts=counts,pos=pos))
}
