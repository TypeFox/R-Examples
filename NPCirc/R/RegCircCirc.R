RegCircCirc<-function(x,y,t,bw,method){
	n <- length(x)
	xt<-outer(x,t,"-")
	if (method=="NW"){
		weights <- 1/(2*pi*besselI(bw,0))*exp(bw*cos(xt))
	}else if (method=="LL"){
		weights <- 1/(2*pi*besselI(bw,0))*exp(bw*cos(xt))
		sinxt<-sin(xt)
		weights<-t(t(n^(-1)*weights)*(colSums(weights*sinxt^2)-t(sinxt)*colSums(weights*sinxt)))
	}
	g1<-colMeans(sin(y)*weights)
	g2<-colMeans(cos(y)*weights)
	fhat <- atan2(g1,g2)
	return(fhat)
}