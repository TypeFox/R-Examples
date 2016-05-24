LikConstantn <- function(lambda,mu,sampling,x,root=1){
		l<-lambda
		m<-mu
		x<-sort(x,decreasing=TRUE)
		t<-x
		rho<-sampling
		lik<- 0+root*(-log(length(t)))
		for (i in 2:length(t)){
		lik<- lik+log(l*p1(t[i],l,m,rho))-log(qhelp(t[1],l,m,rho))}
	 -lik
}