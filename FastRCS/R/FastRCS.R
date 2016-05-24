FRCSnumStarts<-function(p,gamma=0.99,eps=0.5){
	if(p>25)	stop("p too large.")
	if(gamma>=1)	stop("gamma should be smaller than 1.")
	if(gamma<=0)	stop("gamma should be larger than 0.")
	if(eps>0.5)	stop("eps should be smaller than 1/2.")
	if(eps<=0)	stop("eps should be larger than 0.")	
	ns0<-ceiling(log(1-gamma)/log(1-(1-(eps))^(p+1)))
	ord<-10^floor(log10(ns0))
	max(100,ceiling(ns0/ord)*ord)
}
FastRCS<-function(x,y,nSamp=NULL,alpha=0.5,seed=1,intercept=1){#x<-x0;y<-y0;nSamp<-ns;alpha<-0.5
	k1<-25;k0<-25;J<-3;
	m1<-"seed should be an integer in [0,2**31]."
	if(!is.null(seed)){
		if(!is.finite(seed))		stop(m1)
		if(!is.numeric(seed))		stop(m1)
		if(seed<0)			stop(m1)
		if(is.na(as.integer(seed)))	stop(m1)
	}
	seed<-as.integer(seed)+1
	x<-data.matrix(x)
	y<-data.matrix(y)
	na.x<-complete.cases(cbind(x,y))
	if(!is.numeric(alpha))	stop("alpha should be numeric")
	if(alpha<0.5 | alpha>=1)stop("alpha should be in (0.5,1(.")
	if(sum(na.x)!=nrow(x))  stop("Your data contains NA.")
	if(nrow(x)<(5*ncol(x))) stop("n<5p. You need more observations")
	if(intercept){
		cx<-cbind(1,x)
	} else {
		cx<-x
	}
	n<-nrow(cx)
	if(nrow(unique(cbind(x,y)))<n)	stop("Your dataset contains duplicated rows. Please remove them.") 
	p<-ncol(cx)
	if(p<2)			stop("Univariate RCS is not implemented.")
	if(p>25)		stop("FastRCS only works for dimensions <=25.")
	if(is.null(nSamp)) 	nSamp<-FRCSnumStarts(p,eps=(1-alpha)) 
	h1<-min(n-1,quanf(n=n,p=p,alpha=alpha))
	h0<-quanf(n=n,p=p,alpha=0.5);
	Dp<-rep(1.00,n);
	k0<-max(k0,p+2);
	k1<-max(k1,p+2);
	objfunC<-1e3;
	n1<-rep(0,h0);
	n2<-rep(0,h1)
	icandid<-1:n-1
	ni<-length(icandid)
	fitd<-.C("fastrcs",
		as.integer(nrow(cx)),	#1
		as.integer(ncol(cx)),	#2
		as.integer(k0),		#3
		as.single(cx),		#4
		as.single(y),		#5
		as.integer(k1),		#6
		as.single(Dp),		#7
		as.integer(nSamp),	#8
		as.integer(J),		#9
		as.single(objfunC),	#10
		as.integer(seed),	#11
		as.integer(icandid),	#12
		as.integer(ni),		#13
		as.integer(n1),		#14
		as.integer(n2),		#15
		as.integer(h0),		#16
		as.integer(h1),		#17
		PACKAGE="FastRCS")
	outd<-as.numeric(fitd[[7]])
	if(is.nan(outd)[1])	stop("too many singular subsets encoutered!")
	best<-fitd[[15]]#which(outd<=median(outd))
	weit<-as.numeric((1:n)%in%best)
	rawC<-lm(y~x,weights=weit)	
	weit<-as.numeric(abs(rawC$resid)/quantile(abs(rawC$resid),h1/n)/qnorm(1-alpha/2)<=sqrt(qchisq(0.975,df=1)))
       	FinalFit<-lm(y~x,weights=weit);
	FinalFit$scale<-sqrt(sum(FinalFit$resid[weit==1]**2)/FinalFit$df)
	A1<-list(alpha=alpha,nSamp=nSamp,obj=as.numeric(fitd[[10]]),rawBest=fitd[[14]],rawDist=fitd[[7]],  
	best=which(weit==1),coefficients=FinalFit$coef,fitted.values=FinalFit$fitted.values,residuals=FinalFit$residuals,  rank=FinalFit$rank,weights=FinalFit$weights,df.residual=FinalFit$df.residual,scale=FinalFit$scale)
	  class(A1) <-"FastRCS"
	  return(A1)
}
quanf<-function(n,p,alpha)	return(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alpha))
plot.FastRCS<-function(x,col="black",pch=16,...){
  plot(abs(x$resid)/x$scale,col=col,pch=pch,ylab="Robust standardized residuals",xlab="Index")
  abline(h=sqrt(qchisq(0.975,df=1)),col="red",lty=2)
}
