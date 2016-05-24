FPCSnumStarts<-function(p,gamma=0.99,eps=0.5){
	if(p>25)	stop("p too large.")
	if(gamma>=1)	stop("gamma should be smaller than 1.")
	if(gamma<=0)	stop("gamma should be larger than 0.")
	if(eps>0.5)	stop("eps should be smaller than 1/2.")
	if(eps<=0)	stop("eps should be larger than 0.")	
	ns0<-ceiling(log(1-gamma)/log(1-(1-(eps))^(p+1)))
	ord<-10^floor(log10(ns0))
	max(100,ceiling(ns0/ord)*ord)
}
FastPCS<-function(x,nSamp=NULL,alpha=0.5,seed=1){
	k1<-25;k0<-25;J<-3;
	m1<-"seed should be in [0,2**31]."
	if(!is.null(seed)){
		if(!is.finite(seed))		stop(m1)
		if(!is.numeric(seed))		stop(m1)
		if(seed<0)			stop(m1)
		if(is.na(as.integer(seed)))	stop(m1)
	}
	seed<-as.integer(seed)
	x<-data.matrix(x)
	na.x<-complete.cases(x)
	if(!is.numeric(alpha))	stop("alpha should be numeric")
	if(alpha<0.5 | alpha>=1)stop("alpha should be in (0.5,1(.")
	if(sum(na.x)!=nrow(x))  stop("Your data contains NA.")
	if(nrow(x)<(5*ncol(x))) stop("n<5p. You need more observations")
	n<-nrow(x)
	p<-ncol(x)
	if(p<2)		stop("Univariate PCS is not implemented.")
	if(p>25)		stop("FastPCS only works for dimensions<=25.")
	if(is.null(nSamp)) 	nSamp<-FPCSnumStarts(p,eps=(1-alpha)) 
	h1<-min(quanf(n=n,p=p,alpha=alpha),n-1)
	h0<-quanf(n=n,p=p,alpha=0.5)
	Dp<-rep(1.00,n);
	k0<-max(k0,p+1);
	k1<-max(k1,p+1);
	objfunC<-1e3;
	n1<-rep(0,h0);
	n2<-rep(0,h1)
	icandid<-1:n-1
	ni<-length(icandid)
	fitd<-.C("fastpcs",
		as.integer(nrow(x)),	#1
		as.integer(ncol(x)),	#2
		as.integer(k0),		#3
		as.single(x),		#4
		as.integer(k1),		#5
		as.single(Dp),		#6
		as.integer(nSamp),	#7
		as.integer(J),		#8
		as.single(objfunC),	#9
		as.integer(seed),	#10
		as.integer(icandid),	#11
		as.integer(ni),		#12
		as.integer(n1),		#13
		as.integer(n2),		#14
		as.integer(h0),		#15
		as.integer(h1),		#16
		PACKAGE="FastPCS")
	outd<-as.numeric(fitd[[6]])
	if(is.nan(outd)[1])	stop("too many singular subsets encoutered!")	
	best<-as.numeric(fitd[[13]])
	invc<-svd(cov(x[best,]),nu=0,nv=0)
	if(min(invc$d)>1e-8){
		best<-as.numeric(fitd[[14]])
		solV<-chol2inv(chol(cov(x[best,])));
		stds<-mahalanobis(x,colMeans(x[best,]),solV,inverted=TRUE)
		thr0<-qchisq(0.975,df=p)/qchisq(h1/n,df=p)*quantile(stds,h1/n)
		weit<-(stds<=thr0)
		solV<-chol2inv(chol(cov(x[weit,])));
		stds<-sqrt(mahalanobis(x,colMeans(x[weit,]),solV,inverted=TRUE))
	} else {
		"%ni%"<-Negate("%in%") 
		stds<-as.numeric((1:n)%ni%best)
		weit<-!stds
		print("FastPCS has found n/2 observations on a subspace.")
	}
	A1<-list(alpha=alpha,nSamp=nSamp,obj=as.numeric(fitd[[9]]),rawBest=as.numeric(fitd[[13]]),weight=weit,best=best,center=colMeans(x[best,]),cov=cov(x[best,]),distance=stds)
  	class(A1)<-"FastPCS"
	return(A1)
}
quanf<-function(n,p,alpha)	return(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alpha))
plot.FastPCS<-function(x,col="black",pch=16,...){
  plot(x$distance,col=col,pch=pch,ylab="Robust statistical distance",xlab="Index")
  abline(h=sqrt(qchisq(0.975,df=length(x$center))),col="red",lty=2)
}
