
p0est<-function(pvalue,lambda=seq(0,0.90,0.05), pi0.method="smoother",smooth.df=3, smooth.log.pi0=FALSE,fdr.level=NULL){
	z<-k.checkings.pvalue(pvalue=pvalue,p0=NULL)
	zx<-try(qvalue(p=z$pvalue, lambda=lambda, pi0.method=pi0.method, fdr.level=fdr.level,
		smooth.df=smooth.df, smooth.log.pi0=smooth.log.pi0))
	if(is_error(zx)){1}
	else{zx$pi0}
}


##---------------
k.checkings.pvalue<-function(pvalue,p0=NULL){
	pvalue<-MakeNames(pvalue)
	stopifnot(is.null(p0)||is_prob(p0))
	orig.pvalue<-pvalue
	pvalue<-pvalue[is.finite(pvalue)]
	stopifnot(all(is_prob(pvalue)))
	list(pvalue=pvalue,orig.pvalue=orig.pvalue,p0=p0)
}
#---
k.bbe<-function(pvalue,p0=NULL,robust=FALSE,monotonic=FALSE,...){x<-pvalue;method<-'BBE'
	
	if(is.null(p0)){p0<-try(p0est(pvalue=x,...))} 
	if(is_error(p0)){p0<-1;message('p0 estimation failed, setting p0=1');method<-'BBE1'}
	stopifnot(is_prob(p0))
	if (robust==FALSE)
		{z<-try(estimated.LFDR (object=as(x,'Numeric'),achieved.BFDR.fun =b.notrobust,P0=p0,monotonic=monotonic,verbose=FALSE))}
	else if (robust==TRUE)
		{z<-try(estimated.LFDR (object=as(x,'Numeric'),achieved.BFDR.fun =b.robust,P0=p0,monotonic=monotonic,verbose=FALSE))}z<-as(z,"numeric")
	list(lfdr=z,p0=p0,method=method)
}



#-----
k.rvalue<-function(pvalue,robust = FALSE,...){
	zo<-k.qvalue(pvalue=pvalue,robust = robust,...)
	qval<-zo$lfdr;p0<-zo$p0
	if(length(qval)==0){zo<-list(lfdr=NULL,p0=NULL,qvalue=NULL);return(zo)}

	m<-length(qval)
	u <- order(qval)
	rr<-qval[u]
	k.fun<-function(i){
		if (2*i<=m){rr[i]<-min(qval[u[2*i]],1)}
		else{rr[i]<-1}
	rr[i]
	}
	nr<-vapply(1:m,FUN=k.fun,FUN.VALUE=numeric(1))
	names(nr)<-names(rr)
	list(lfdr=nr,p0=p0,qvalue=qval)	
}
k.qvalue<-function(pvalue,robust = FALSE,...){
	qvalx<-try(qvalue(p=pvalue,robust = robust, ...))
	
	if(is_error(qvalx)|| is(qvalx,'numeric')){qx<-list(lfdr=NULL,p0=NULL);return(qx)}
	else{
		assert.is(qvalx,"qvalue")
		p0<-qvalx$pi0
		qval<-qvalx$qvalues}

	
	list(lfdr=qval,p0=p0)
}

#---------

PHATs.pvalue<-function(lfdr.fun='rvalue',pvalue,p0=NULL,robust=FALSE,monotonic=FALSE,...){
	assert.is(pvalue,'numeric')
		
	z<-k.checkings.pvalue(pvalue=pvalue,p0=p0)
	#---------------
	if(lfdr.fun%in%c('nqvalue','qvalue','q')){
		method<-'QVALUES'
		zo<-k.qvalue(pvalue=z$pvalue,robust=robust,...)
		infox<-c(list(...),list(robust=robust))
		infox$method<-method
		
	}
	else if(lfdr.fun%in%c('rvalue','r')){
		method<-'RVALUES'
		zo<-k.rvalue(pvalue=z$pvalue,robust=robust,...)
		infox<-c(list(...),list(robust=robust))
		infox$method<-method
	}
	else if(lfdr.fun%in%c('lfdr.bbe','bbe')){
		zo<-k.bbe(pvalue=z$pvalue,p0=z$p0,robust=robust,monotonic=monotonic,...)
		infox<-c(list(...),list(p0=p0,robust=robust,monotonic=monotonic))
		infox$method<-method<-zo$method
	}
	else if(lfdr.fun%in%c('lfdr.bbe1','bbe1')){z$p0<-1
		zo<-k.bbe(pvalue=z$pvalue,p0=z$p0,robust=robust,monotonic=monotonic,...)
		infox<-c(list(...),list(p0=p0,robust=robust,monotonic=monotonic))
		infox$method<-method<-'BBE1'
	}
	#--------------------
	zo<-new_est.lfdr.pvalue(LFDR.hat=zo$lfdr,p0.hat=zo$p0,pvalue=z$orig.pvalue,method=method,info=c(uniquex(c(zo$info,infox)))	)
	est2list(zo)
	
}


#
rvalue<-function(pvalue,robust=FALSE,...){
	PHATs.pvalue(lfdr.fun='rvalue',pvalue=pvalue,p0=NULL,robust=robust,monotonic=FALSE,...)	
}
nqvalue<-function(pvalue,robust=FALSE,...){
	PHATs.pvalue(lfdr.fun='nqvalue',pvalue=pvalue,p0=NULL,robust=robust,monotonic=FALSE,...)	
}
#
lfdr.bbe1<-function(pvalue,robust=FALSE,monotonic=FALSE,...){
	PHATs.pvalue(lfdr.fun='lfdr.bbe1',pvalue=pvalue,p0=1,robust=robust,monotonic=monotonic,...)	
}
lfdr.bbe<-function(pvalue,p0=NULL,robust=FALSE,monotonic=FALSE,...){
	PHATs.pvalue(lfdr.fun='lfdr.bbe',pvalue=pvalue,p0=p0,robust=robust,monotonic=monotonic,...)	
}
#----
