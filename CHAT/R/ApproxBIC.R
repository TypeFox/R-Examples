ApproxBIC <- function(Y,fit.x,fit.y,n){
	L<-0
	for(y in Y){
		vv1<-which(fit.x<=y)
		vv2<-which(fit.x>=y)
		fhat<-mean(fit.y[vv1[length(vv1)]],fit.y[vv2[1]])
		L<-L-2*log(fhat)
	}
	return(L+n*log(length(Y)))
}
