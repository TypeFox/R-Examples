cyclic.covering<-function(mat,vec.primes){
	ia<-as.data.frame(ia.samp(ncol(mat)))
	not.covered<-rowSums(as.matrix(ia)%*%t(mat)==0)
	ia<-ia[not.covered==0,]
	rowS<-rowSums(ia)
	min.rowS<-which(rowS==min(rowS))
	ia<-ia[min.rowS,]
	list.primes<-vector("list",nrow(ia))
	for(i in 1:nrow(ia))
		list.primes[[i]]<-c(vec.primes,colnames(mat)[ia[i,]==1])
	list.primes
}

