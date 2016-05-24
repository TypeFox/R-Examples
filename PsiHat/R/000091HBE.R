

#============for LFDRs==========

k.valid.x<-function(object,x="pvalue"){
	dx<-object@LFDR.hat
	p0<-object@p0.hat
	stx<-eval(parse(text=paste("object@",x,sep="")))
	ix<-object@info
	nvar<-length(dx)
	stopifnot(length(names(dx))>0)
	stopifnot(length(p0)%in%c(1,nvar))
	stopifnot(length(stx)%in%c(0,nvar))
	if(length(stx)==nvar){stopifnot(identical(names(stx),names(dx)))}
	if(!are_prob(list(dx[is.finite(dx)],p0[is.finite(p0)]))){stop('error: LFDR.hat or p0.hat out of [0,1]')}
	if(x%in%c("pvalue")&&!is_prob(stx[is.finite(stx)])){stop('error: pvalue out of [0,1]')}
	TRUE}


setClass("est.lfdr.pvalue", representation(LFDR.hat ="numeric",p0.hat = "numeric",pvalue="numeric",info="list"))
setValidity("est.lfdr.pvalue", function(object) {k.valid.x(object,x="pvalue")})
setClass("est.lfdr.stat", representation(LFDR.hat ="numeric",p0.hat = "numeric", stat="numeric",info="list"))
setValidity("est.lfdr.stat", function(object) {k.valid.x(object,x="stat")})
setClassUnion("est.lfdr.x",c("est.lfdr.pvalue","est.lfdr.stat"))

k.new_est.lfdr.x <- function(LFDR.hat,p0.hat,x=NULL,method=NULL,xtype=c("pvalue","stat"),info=list()){
	
	err<-FALSE;xtype<-xtype[1]
	if(!is(info,'list')){info<-list()}
	assert.is(LFDR.hat,c('NULL',"numeric","logical",'try-error'))
	assert.is(x,c("numeric","logical"))
	stopifnot(length(x)>=length(LFDR.hat))
	stopifnot(xtype%in%c("pvalue","stat"))
	if(is_error(LFDR.hat)){err<-TRUE
		message(method," estimation failed ",Sys.time())
		LFDR.hat<-rep(as.numeric(NA),length(x))
		names(LFDR.hat)<-names(x)
		method<-paste('FAILED-',method,sep="")
		p0.hat<-as.numeric(NA)}
	
	nvar<-length(x)
	names(p0.hat)<-NULL
	
	
	if(is(x,"logical")){xx<-as.numeric(x);names(xx)<-names(x);x<-xx}
	if(is(LFDR.hat,"logical")){xx<-as.numeric(LFDR.hat);names(xx)<-names(LFDR.hat);LFDR.hat<-xx}
	if(is(p0.hat,"logical")){xx<-as.numeric(p0.hat);names(xx)<-names(p0.hat);p0.hat<-xx}	
	
	assert.are(list(LFDR.hat,p0.hat,x),'numeric')
	stopifnot(length(p0.hat)%in%c(1,nvar))
	x<-MakeNames(x,nmvar="X")
	LFDR.hat<-MakeNames(LFDR.hat,nmvar="X")
	LFDR.hat<-sameAsX_names(y=LFDR.hat,x=x)

	nvar<-length(LFDR.hat)
	stopifnot(length(x)%in%c(nvar))
	stopifnot(identical(names(x),names(LFDR.hat)))

	if(!err){
		p0.hat<-pmin(p0.hat,rep(1,length(p0.hat)))
		stopifnot(all(is_prob(LFDR.hat[is.finite(LFDR.hat)])) &&  all(is_prob(p0.hat[is.finite(p0.hat)])))
	}

	if(is.null(method)){method<-'?'}
	info<-info[!names(info)%in%c('x','y','W','stat','p.value','pvalue')]

	infox<-uniquex(c(list(method=method),info))
	
	if(xtype%in%c("pvalue")){stopifnot(all(is_prob(x[is.finite(x)])))}

	zo<-eval(parse(text=paste('new("est.lfdr.',xtype,'",LFDR.hat=LFDR.hat, p0.hat=p0.hat,',xtype,'=x,info=infox)',sep="")))
	validObject(zo)
	zo
}

new_est.lfdr.stat <- function(LFDR.hat,p0.hat,stat,method=NULL,info=list()){
	if(is_nothing(method)){method<-"HBE?"}
	k.new_est.lfdr.x(LFDR.hat=LFDR.hat,p0.hat=p0.hat,x=stat,method=method,xtype=c("stat"),info=info)}
new_est.lfdr.pvalue <- function(LFDR.hat,p0.hat,pvalue,method=NULL,info=list()){
	if(is_nothing(method)){method<-"?"}
	k.new_est.lfdr.x(LFDR.hat=LFDR.hat,p0.hat=p0.hat,x=pvalue,method=method,xtype=c("pvalue"),info=info)}



setMethod("[", signature(x = "est.lfdr.x", i = "ANY", j = "missing"), function(x, i, j, drop){
	z<-x;xtype<-slotNames(x);xtype<-xtype[xtype%in%c("pvalue","stat")][1]
	nvar<-length(x@LFDR.hat)
	z@LFDR.hat <- x@LFDR.hat[i];names(z@LFDR.hat)<-names(x@LFDR.hat)[i]
	eval(parse(text=paste("z@",xtype," <- x@",xtype,"[i]",sep="")))
	eval(parse(text=paste("names(z@",xtype,")<-names(x@LFDR.hat)[i]",sep="")))
	stopifnot(validObject(z))
	z})

setReplaceMethod("names", signature(x="est.lfdr.x",value="character"),function(x, value){
  nvar<-length(x@LFDR.hat)
  xtype<-slotNames(x);xtype<-xtype[xtype%in%c("pvalue","stat")][1]
  stopifnot(length(value)==nvar)
  names(x@LFDR.hat)<-value
  eval(parse(text=paste("names(x@",xtype,") <- value",sep="")))
  stopifnot(validObject(x))
	x})
                 
setMethod("names", signature(x = "est.lfdr.x"), function(x){names(x@LFDR.hat)})
setMethod("as.numeric", signature(x="est.lfdr.x"),function(x) {x@LFDR.hat})
setMethod("length", signature(x = "est.lfdr.x"), function(x){length(as.numeric(x))})


setMethod("is_prob",signature(object = "est.lfdr.stat" ), function(object,tolerance = 1e-3, ...){
         is_prob(object@LFDR.hat,tolerance =tolerance,...) && is_prob(object@p0.hat,tolerance =tolerance,...)})

setMethod("is_prob",signature(object = "est.lfdr.pvalue" ), function(object,tolerance = 1e-3, ...){
         is_prob(object@pvalue,tolerance =tolerance,...) &&is_prob(object@LFDR.hat,tolerance =tolerance,...) && is_prob(object@p0.hat,tolerance =tolerance,...)})
##------------------------------


mylocfdr <- locfdr

#----------------------

k.locfdr2est.lfdr.pvalue<-function(x,xo,nulltype,info=list()){
	if(is_error(x)){return(list())}
	assert.is(x,c('list','locfdr.list'))
	assert.is(xo,c('numeric'))
	stopifnot(length(xo)==length(x$fdr))
	names(x$fdr)<-names(xo)
rname <- locfdr.rname(nulltype = nulltype)
p0 <- x$fp0[rname, "p0"]
if(nulltype%in%0){method<-'HBEA'}
else if(nulltype%in%1){method<-'HBEE'}
else if(nulltype%in%c(2,3)){method<-paste('HBEE',nulltype,sep='')}
else{stop('bad nulltype')}
	list(lfdr=x$fdr,p0=p0,method=method,info=uniquex(c(list(method=method),info)))
}
k.hbe.method<-function(nulltype){
	if(nulltype%in%0){method<-'HBEA'}
else if(nulltype%in%1){method<-'HBEE'}
else if(nulltype%in%c(2,3)){method<-paste('HBEE',nulltype,sep='')}
else{stop('bad nulltype')}
method
}
#----------------
k.checkings.hbe<-function(pvalue=NULL,stat=NULL,bre=120, df=7, nulltype=0,plot=0){
	assert.are(list(bre,df,nulltype,plot),'numeric')
	stopifnot(nulltype%in%c(0:3))
	stopifnot(plot%in%c(0:4))
	stopifnot(length(stat[is.finite(stat)])>0|length(pvalue[is.finite(pvalue)])>0)
       	
        if(length(pvalue)>0 && length(stat)>0){message("HBE from statistics vector");pvalue<-NULL}

	orig.stat<-orig.pvalue<-NULL
        
	if (length(stat) > 0){orig.stat<-stat <- MakeNames(stat)}
	if (length(pvalue[is.finite(pvalue)]) > 0){
		stopifnot(all(is_prob(pvalue[is.finite(pvalue)])))
		pvalue <- MakeNames(pvalue);orig.pvalue<-pvalue
		orig.stat<-stat<-qnorm(pvalue)}
	
		
	stopifnot(any(is.finite(stat))&&any(is.finite(orig.stat)))
	stopifnot(identical(names(stat),names(orig.stat)))

	
	ind<-which(is.finite(stat))
	stat<-stat[ind]
	
		
	method<-k.hbe.method(nulltype=nulltype)
	list(pvalue=pvalue,stat=stat,orig.pvalue=orig.pvalue,orig.stat=orig.stat,bre=bre,df=df,nulltype=nulltype,plot=plot,method=method)
}
#-

k.empNull<-function(stat,bre, df, nulltype, plot=0,...){
	#x<-stat;arg.cv<-list(...);arg.cv<-arg.cv[names(arg.cv)%in%c('s3FUN', 'arglis')]
	#cv<-do.call(pval2cval,c(list(pvalue=pvalue),arg.cv))   ; cv<-x
	
	z<-try(nempiricalNull(object=stat, plot=plot, nsilence = 0, silent = NULL,
			     call.browser = FALSE, cvalue.arglis = NULL, verbose = TRUE, max.p0 = 1,
			     bre=bre, df=df, nulltype=nulltype,...))
	if(is_error(z)){zo<-list()}
	else{zo<-list(lfdr=z@s3$fdr,p0.hat=as(z@p0,'numeric'),empNull=z)}
zo
}

k.elfdr<-function(stat,bre, df, nulltype, plot=0,...){x<-stat
	zx<-try(k.empNull(stat=x,bre=bre, df=df, nulltype=nulltype, plot=plot,...),silent=TRUE)
	if(is_error(zx)){zo<-list();return(zo)}
	plots<-ifelse(plot==0,FALSE,TRUE)
	z<-try(expected.lfdr(object=zx$empNull, call.plot =  plots))
	if(is_error(z)){zo<-list()}
	else{zo<-list(lfdr=CorrectLimits(z),p0.hat=as(zx$ empNull@p0,'numeric'))}
	zo}
#----------------
PHATs.stat<-function(lfdr.fun='lfdr.hbea',stat=NULL,pvalue=NULL,plot=0,nulltype=1,bre=120, df=7,...){
	
	
	z<-k.checkings.hbe(pvalue=pvalue,stat=stat,bre=bre, df=df, nulltype=nulltype,plot=plot)

	#--------------------
	if(lfdr.fun%in%c('lfdr.hbea','hbea')){
		z$nulltype<-0
		zx<-try(locfdr(zz=z$stat, bre=z$bre, df=z$df, nulltype=z$nulltype,  plot=z$plot,...))
		zo<-k.locfdr2est.lfdr.pvalue(x=zx,xo=z$stat,nulltype=z$nulltype,info=list())
                method<-z$method#"HBEA"
	}
	else if(lfdr.fun%in%c('lfdr.assumedNull','lfdr.assNull','assumedNull','assNull','asNull')){
		z$nulltype<-0
		zo<-k.empNull(stat=z$stat,bre=z$bre, df=z$df, nulltype=z$nulltype, plot=z$plot,...)
                method<-"AssumedNull"
	}
	else if(lfdr.fun%in%c('lfdr.hbee','hbee')){#z$nulltype<-1
		zx<-try(locfdr(zz=z$stat, bre=z$bre, df=z$df, nulltype=z$nulltype, plot=z$plot,...))
		zo<-k.locfdr2est.lfdr.pvalue(x=zx,xo=z$stat,nulltype=z$nulltype,info=list())
		method<-z$method
	}
	else if(lfdr.fun%in%c('lfdr.empiricalNull','lfdr.empNull','empNull')){#z$nulltype<-1
		zo<-k.empNull(stat=z$stat,bre=z$bre, df=z$df, nulltype=z$nulltype, plot=z$plot,...)
		method<-gsub(pattern='HBEE', replacement='EmpiricalNull', x=z$method, ignore.case = FALSE, fixed=TRUE)
	}
	else if(lfdr.fun%in%c('lfdr.elfdr')){#z$nulltype<-1
		zo<-k.elfdr(stat=z$stat,bre=z$bre, df=z$df, nulltype=z$nulltype, plot=z$plot,...)
		method<-'ELFDR';
		
	}
	else if(lfdr.fun%in%c('lfdr.hbe','hbe')){
		method<-z$method
		zx<-try(locfdr(zz=z$stat, bre=z$bre, df=z$df, nulltype=z$nulltype,  plot=z$plot,...))
		zo<-k.locfdr2est.lfdr.pvalue(x=zx,xo=z$stat,nulltype=z$nulltype,info=list())
		method<-"locfdr";
	}
	else{stop('bad lfdr.fun')}
	
	infox<-c(list(...),z[c("nulltype","bre","df")])
	infox$method<-method
	#--------------------
	zo<-new_est.lfdr.stat(LFDR.hat=zo$lfdr,p0.hat=zo$p0,stat=z$orig.stat,method=method,info=c(uniquex(c(zo$info,infox))))
	est2list(zo)
}

#-------
lfdr.elfdr<-function(stat= NULL, pvalue = NULL, nulltype = 1, bre = 120, df = 7, plot = 0, ...){
		     
	PHATs.stat(lfdr.fun='lfdr.elfdr',pvalue=pvalue, stat=stat,plot=plot,nulltype=nulltype,bre=bre, df=df,...)
	
}

lfdr.empiricalNull<-function(stat= NULL, pvalue = NULL,  nulltype = 1, bre = 120, df = 7, plot = 0, ...){

	PHATs.stat(lfdr.fun='lfdr.empiricalNull',pvalue=pvalue, stat=stat,plot=plot,nulltype=nulltype,bre=bre, df=df,...)
	
}
lfdr.assumedNull<-function(stat= NULL, pvalue = NULL,  bre = 120, df = 7, plot = 0, ...){
	nulltype<-0
	PHATs.stat(lfdr.fun='lfdr.assumedNull', pvalue=pvalue, stat=stat,plot=plot,nulltype=nulltype,bre=bre, df=df,...)
	
}
lfdr.hbee<-function(stat= NULL, pvalue = NULL,  nulltype = 1, bre = 120, df = 7, plot = 0, ...){
	PHATs.stat(lfdr.fun='lfdr.hbee', pvalue=pvalue, stat=stat,plot=plot,nulltype=nulltype,bre=bre, df=df,...)
	
}
lfdr.hbea<-function(stat= NULL, pvalue = NULL,  nulltype = 1, bre = 120, df = 7, plot = 0, ...){
	nulltype<-0
	PHATs.stat(lfdr.fun='lfdr.hbea', pvalue=pvalue, stat=stat,plot=plot,nulltype=nulltype,bre=bre, df=df,...)
	
}
lfdr.hbe<-function(stat= NULL, pvalue = NULL, nulltype = 1, bre = 120, df = 7, plot = 0, ...){
	PHATs.stat(lfdr.fun='lfdr.hbee', pvalue=pvalue, stat=stat,plot=plot,nulltype=nulltype,bre=bre, df=df,...)
	
}
###-------------------------------------
