identifyMonomorphism<-function(x){
	mx<-max(x,na.rm=TRUE)
	checkCatMat(x,mx,matname="x")
	n<-rowSums(!is.na(x))
	rs<-numeric(nrow(x))
	for(i in 1:mx){
		tmp<-rowSums(x==i,na.rm=TRUE)==n
		rs<-rs+tmp
	}
	which(rs>0)
}

