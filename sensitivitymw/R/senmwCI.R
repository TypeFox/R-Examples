senmwCI <- function(y,gamma=1,method=NULL,inner=0,trim=3,lambda=1/2,m1=1,m2=1,m=1,alpha=0.05,one.sided=TRUE,tol=NULL,interval=NULL,detail=FALSE){
	if (is.vector(y)) {
		y <- y[!is.na(y)]
		treat <- y/2
		cont <- (-y/2)
		y <- cbind(treat, cont)
	}
	if (m2<m) warning("Redescending scores, m2<m, may not yield sensible confidence intervals and estimates")
	stopifnot(m2==m)
	if (!is.null(method)) {if (method=="l"){
		warning("method=l with redescending scores may not yield sensible confidence intervals and estimates")
		stopifnot(method!="l")
		}
	}
	stopifnot ((alpha>0)&(alpha<1))
	if (!one.sided) alpha <- alpha/2
	
	funcCI<-function(tau){
		target<-qnorm(1-alpha)
		ntau<-length(tau)
		o<-rep(NA,ntau)
		for (i in 1:ntau){
			dev<-senmw(y,gamma=gamma,method=method,inner=inner,trim=trim,lambda=lambda,m1=m1,m2=m2,m=m,tau=tau[i])$deviate
			o[i]<-(dev-target)
		}
		o
	}
	
	funcCI2<-function(tau){
		target<-qnorm(1-alpha)
		ntau<-length(tau)
		o<-rep(NA,ntau)
		for (i in 1:ntau){
			dev<-senmw(-y,gamma=gamma,method=method,inner=inner,trim=trim,lambda=lambda,m1=m1,m2=m2,m=m,tau=tau[i])$deviate
			o[i]<-(dev-target)
		}
		o
	}
		
	funcEST<-function(tau){
		target<-0
		ntau<-length(tau)
		o<-rep(NA,ntau)
		for (i in 1:ntau){
			dev<-senmw(y,gamma=gamma,method=method,inner=inner,trim=trim,lambda=lambda,m1=m1,m2=m2,m=m,tau=tau[i])$deviate
			o[i]<-(dev-target)
		}
		o
	}
	
	funcEST2<-function(tau){
		target<-0
		ntau<-length(tau)
		o<-rep(NA,ntau)
		for (i in 1:ntau){
			dev<-senmw(-y,gamma=gamma,method=method,inner=inner,trim=trim,lambda=lambda,m1=m1,m2=m2,m=m,tau=tau[i])$deviate
			o[i]<-(dev-target)
		}
		o
	}
	
	tr<-as.vector(y[,1])
	ct<-as.vector(y[,-1])
	mx<-max(tr)-min(ct)
	mn<-min(tr)-max(ct)
	if (is.null(interval)) interval<-c(mn,mx)
	else {
		stopifnot(length(interval)==2)
	}
	stopifnot(interval[2]>interval[1])
	interval2<-c(-interval[2],-interval[1])
	if (is.null(tol)) tol<-((max(interval)-min(interval)))/500000
	else stopifnot(tol>0)
	vCI<-uniroot(funcCI,interval=interval,tol=tol)
	vEST<-uniroot(funcEST,interval=interval,tol=tol)
	if (!one.sided) {vCI2<-uniroot(funcCI2,interval=interval2,tol=tol)}
	vEST2<-uniroot(funcEST2,interval=interval2,tol=tol)
	min.estimate<-vEST$root
	max.estimate<-(-vEST2$root)
	min.lowerCI<-vCI$root
	PointEstimate<-c(min.estimate,max.estimate)
	names(PointEstimate)<-c("minimum","maximum")
	if (one.sided){CI<-c(min.lowerCI,Inf)}
	else {
		max.upperCI<-(-vCI2$root)
		CI<-c(min.lowerCI,max.upperCI)
	}
	names(CI)<-c("minimum","maximum")
	if (detail) {list(PointEstimate=PointEstimate,Confidence.Interval=CI,search.interval=interval,tolerance=tol)}
	else {
		rd<-floor(log10(1/tol))-1
		list(PointEstimate=round(PointEstimate,rd),Confidence.Interval=round(CI,rd))}
}