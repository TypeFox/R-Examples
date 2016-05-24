bray.part<-function(x) {
	x <- as.matrix(x)
	result<-matrix(nrow=nrow(x),ncol=nrow(x))
	rownames(result)<-rownames(x)
	colnames(result)<-rownames(x)

	for(i in 1:nrow(x)) {
    	for(j in i:nrow(x)) {
		A<-sum(pmin(x[i,],x[j,]))
		B<-sum(x[i,])-sum(pmin(x[i,],x[j,]))
		C<-sum(x[j,])-sum(pmin(x[i,],x[j,]))
		result[i,j]<-min(B,C)/(A+min(B,C))
		result[j,i]<-(B+C)/(2*A+B+C)
		}
		}
	bray<-as.dist(result)
	bray.bal<-as.dist(t(upper.tri(result)*result))
	bray.gra=bray-bray.bal

	results<-list( bray.bal=bray.bal, bray.gra=bray.gra, bray=bray)
	return(results)
}