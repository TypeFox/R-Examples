LikConstant <- function(lambda,mu,sampling,x,root=1,survival=1){
		l<-lambda
		m<-mu
		x<-sort(x,decreasing=TRUE)
		t<-x
		rho<-sampling
		lik<- (root+1)*log(p1(t[1],l,m,rho))
		for (i in 2:length(t)){
		lik<- lik+log(l*p1(t[i],l,m,rho))}
		if (survival == 1){lik<-lik-(root+1)*log(1-p0(t[1],l,m,rho))}
	-lik
}