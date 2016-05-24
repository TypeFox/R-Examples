check4Monomorphism<-function(x){
	mx<-max(x,na.rm=TRUE)
	n<-rowSums(!is.na(x))
	for(i in 1:mx){
		rs<-rowSums(x==i,na.rm=TRUE)==n
		if(any(rs))
			stop("All variables must show at least two levels.")
	}
}


	