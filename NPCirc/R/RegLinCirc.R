RegLinCirc<-function(x,y,t,bw,method){
	n <- length(x)
	xt<-outer(x,t,"-")
	if (method=="NW"){
		weights <- dnorm(xt,0,bw)
	}else if (method=="LL"){
		weights <- dnorm(xt,0,bw)
		weights <- t(t(n^(-1)*weights)*(colSums(weights*xt^2) -t(xt)*colSums(weights*xt,2)))
	}
	g1<-colMeans(sin(y)*weights)
	g2<-colMeans(cos(y)*weights)
	fhat <- atan2(g1,g2)
	return(fhat)
}