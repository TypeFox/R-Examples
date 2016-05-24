FHCSnumStarts<-function(q,gamma=0.99,eps=0.5){
	if(q>25)	stop("q too large.")
	if(gamma>=1)	stop("gamma should be smaller than 1.")
	if(gamma<=0)	stop("gamma should be larger than 0.")
	if(eps>0.5)	stop("eps should be smaller than 1/2.")
	if(eps<=0)	stop("eps should be larger than 0.")	
	ns0<-ceiling(log(1-gamma)/log(1-(1-(eps))^(q+1)))
	ord<-10^floor(log10(ns0))
	ceiling(ns0/ord)*ord
}
FastHCS<-function(x,nSamp=NULL,alpha=0.5,q=10,seed=1){#x<-x0;nSamp<-100;alpha<-0.5;q=5;seed=NULL
	q0<-q1<-25;J<-5;
	m1<-"seed should be an integer in [0,2**31]."
	if(!is.null(seed)){
		if(!is.finite(seed))		stop(m1)
		if(!is.numeric(seed))		stop(m1)
		if(seed<0)			stop(m1)
		if(is.na(as.integer(seed)))	stop(m1)
	}
	seed<-as.integer(seed)+1
	x<-data.matrix(x)
	na.x<-complete.cases(x)
	if(!is.numeric(alpha))		stop("alpha should be numeric")
	if(alpha<0.5 | alpha>=1)	stop("alpha should be in (0.5,1(.")
	if(sum(na.x)!=nrow(x))		stop("Your data contains NA.")
	if(ncol(x)>=nrow(x)){
		z0<-FHCSkernelEVD(x,q=ncol(x))
	} else {
		z0<-list()
		class(z0)<-"FHCSkernelEVD"
		z0$rank<-ncol(x)
		z0$center<-colMeans(x)
		z0$scores<-sweep(x,2,z0$center)
	}
	if(z0$rank<2)	stop("The true dimensionality of your data is 1. Univariate FHCS is not implemented.")
	if(z0$rank<q){
		print(paste0("The true dimensionality of your data is ",max(z0$rank)-1))
		q<-max(z0$rank)-1
	}
	n<-nrow(z0$scores)
	p<-ncol(z0$scores)
	if(!is.numeric(q))	stop("q should satisfy q>0.")
	if(q<2)			stop("Univariate FastHCS is not implemented.")
	if(p<4)			stop("For small values of p, use FastPCS.")
	qf<-q
	if(q==2){
				qf<-2
				q<-3
	}
	if(q>=p)		stop("q should satisfy q<p.")
	if(q>=n)		stop("q should satisfy q<n.")
	if(q>25)		stop("FastHCS only works for q<=25.")
	if(is.null(nSamp)) 	nSamp<-FHCSnumStarts(q,eps=(1-alpha)) 
	h0<-quanf(n=n,p=q,alpha=0.5)
	h1<-min(n-1,quanf(n=n,p=q,alpha=alpha))
	Dp<-rep(1.00,n);
	q0<-max(q0,q);
	q1<-max(q1,q);
	objfunC<-1e3;
	icandid<-1:n-1
	ni<-length(icandid)
	rraw<-rep(0,h0)
	rrew<-rep(0,h1)
	qs<-rep(0,q+1)
	sd.d<-rep(0,n)
	fitd<-.C("FastHCS",	
		as.integer(nrow(z0$scores)),	#01
		as.integer(ncol(z0$scores)),	#02
		as.integer(q0),			#03
		as.single(z0$scores),		#04
		as.integer(q1),			#05
		as.integer(q),			#06
		as.integer(nSamp),		#07
		as.integer(J),			#08
		as.single(objfunC),		#09
		as.integer(seed),		#10
		as.integer(icandid),		#11
		as.integer(ni),			#12
		as.integer(rraw),		#13
		as.integer(h0),			#14
		as.single(sd.d),		#15
		as.integer(h1),			#16
		as.integer(qs),			#17
	PACKAGE="FastHCS")
	eStep<-compPcaParams(x=x,fitd=fitd,q=min(qf,fitd[[6]]),z0=z0,seed=seed);#
	A1<-list(rawBest=fitd[[13]],obj=as.numeric(fitd[[9]]),rawDist=fitd[[15]],best=eStep$best,
	loadings=eStep$loadings,eigenvalues=eStep$eigenvalues,center=eStep$center,od=eStep$od,
	sd=eStep$sd,cutoff.od=eStep$cutoff.od,cutoff.sd=eStep$cutoff.sd,scores=eStep$scores,
	qstar=fitd[[17]],flag=eStep$flag,best_alt=eStep$alt,z0=z0)
	class(A1)<-"FastHCS"
	return(A1)
}
compPcaParams<-function(x,fitd,q=NULL,z0=NULL,seed=1){
	n<-nrow(x);p<-ncol(x)
	best1<-which(fitd[[15]]<quantile(fitd[[15]],fitd[[16]]/n));
	est0<-FHCSpsdo(z0=z0,h=length(best1),q=q,seed=seed);
	est0$scores<-sweep(x,2,est0$center)%*%est0$loadings
	est1<-FHCSkernelEVD(x=x,best=best1,q=q)							
	cand<-intersect(est1$best,est0$best)
	kand<-setdiff(est0$best,cand)
	cond0<-mean(log(colMeans(est1$scores[est1$best,]**2)/colVars(est1$scores[cand,])))
	cond1<-max(log(colMeans(est0$scores[cand,]**2)/colVars(est0$scores[kand,])))
	cond<-cond0<cond1
	if(max(colVars(est0$scores[kand,]))<1e-6)	cond<-FALSE
	if(cond){
		flag<-0
		best<-est1$best
		center<-est1$center
		scores<-est1$scores
		loadings<-est1$loadings
		eigenvalues<-est1$eigenvalues
		rank<-est1$rank
		alt<-est0$best
	} else {
		flag<-1
		best<-est0$best
		scores<-est0$scores
		center<-est0$center
		loadings<-est0$loadings
		eigenvalues<-est0$eigenvalues
		rank<-est0$rank
		alt<-est1$best
	}
	dod<-sweep(x,2,center,FUN='-')-scores%*%t(loadings);
	dod<-rowSums(dod*dod)
	dsd<-sweep(scores,2,sqrt(eigenvalues),FUN="/")
	dsd<-rowSums(dsd*dsd)	
	cof<-quantile(dsd,fitd[[16]]/n)**2/qchisq(fitd[[16]]/n,rank)
	eigenvalues<-eigenvalues*cof
	dsd<-dsd/sqrt(cof)
	cod<-sqrt(qnorm(0.975,mean(dod[best]^(2/3)),sd(dod[best]^(2/3))/sqrt(qchisq(length(best)/n,1)))^3)
	csd<-sqrt(qchisq(0.975,rank))
	list(center=center,loadings=loadings,eigenvalues=eigenvalues,flag=flag,
	scores=scores,cutoff.sd=csd,cutoff.od=cod,od=dod,sd=dsd,best=best,alt=alt)
}
quanf<-function(n,p,alpha)	return(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alpha))
plot.FastHCS<-function(x,col="black",pch=16,...){
	SDIND<-x$sd
	ODIND<-x$od
	plot(SDIND,ODIND,col=col,pch=pch,xlab="Robust score distance",ylab="Robust orthogonal distance",main="Robust PCA")
	abline(h=x$cutoff.od,col="red",lty=2)
	abline(v=x$cutoff.sd,col="red",lty=2)
}
FHCSkernelEVD<-function(x,best=NULL,q=NULL){#lightly modified from rrcov:
	p<-ncol(x)
	n<-nrow(x)
	if(is.null(q))		q<-p
	if(is.null(best))	best<-1:n
	center<-colMeans(x[best,])
	x<-sweep(x,2,center,FUN='-',check.margin=FALSE)
	e<-eigen(tcrossprod(x[best,])/(length(best)-1),symmetric=TRUE)
	tolerance<-n*max(e$values)*.Machine$double.eps
	rank<-min(q,sum(e$values>tolerance))
	eigenvalues<-e$values[1:rank]
	loadings<-t((x[best,]/sqrt(length(best)-1)))%*%e$vectors[,1:rank]%*%diag(1/sqrt(eigenvalues))
	loadings<-signFlip(loadings)
	scores<-x%*%loadings	
	A1<-list(loadings=loadings,rank=rank,center=center,scores=scores,eigenvalues=eigenvalues,best=best)
	class(A1)<-"FHCSkernelEVD"
	return(A1)
}#from rrcov:
signFlip<-function(loadings)	apply(loadings,2,function(x) if(x[which.max(abs(x))]<0) -x else x)
FHCSpsdo<-function(z0=NULL,h=NULL,seed=1,q=NULL,ndir=1000){
	if(is(z0,'matrix')){
		z0<-FHCSkernelEVD(z0,q=ncol(z0))
	} else {
		if(is(z0,'FHCSkernelEVD')==FALSE)	stop('z0 should either be a numeric matrix or an object of class FHCSkernelEVD.')
	}
	if(is.null(q))	q<-z0$rank
	n<-nrow(z0$scores)
	p<-ncol(z0$scores)
	if(is.null(h))	h<-ceiling((n+q+1)*0.5)
	fitd<-.C("r_psdo",	
		as.integer(n),			#1
		as.integer(p),			#2
		as.integer(ndir),		#3
		as.single(z0$scores),		#4
		as.single(rep(0,n)),		#5
		as.integer(0),			#6
		as.integer(h),			#7
		as.integer(seed),		#8
		as.integer(1:h),		#9
	PACKAGE="FastHCS")
	Xh.svd<-FHCSkernelEVD(x=z0$scores,best=fitd[[9]],q=q)	
	if(!is.null(z0$loadings)){
		center<-z0$center+Xh.svd$center%*%t(z0$loadings)
	        loadings<-z0$loadings%*%Xh.svd$loadings
	} else {
		center<-z0$center+Xh.svd$center
	        loadings<-Xh.svd$loadings
	}
	eigenvalues<-Xh.svd$eigenvalues
	return(list(loadings=loadings,center=center,eigenvalues=eigenvalues,
	rawDist=fitd[[5]],rank=Xh.svd$rank,best_D=fitd[[9]]))
}
