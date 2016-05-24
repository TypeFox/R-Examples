unMAP <-
function(vec){
	N<-length(vec)
	unmap1<-matrix(0,N,max(vec))
	for(ind in 1:N) unmap1[ind, vec[ind]] <- 1
	unmap1
	}
