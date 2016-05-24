##------------------------------util
print_info<-echo<-FALSE
est2list<-function(x){

	if (is(x,"list")){zx<-lapply(x,est2list);return(zx)}
	slnm<-slotNames(x)

	if(!isS4(x)){return(x)}
 
	z<-lapply(1:length(slnm),function(i){
	     z1<-eval(parse(text=paste("x@",slnm[i],sep="")))	
	     if(!isS4(z1)){return(z1)}
	     else {z1<-est2list(x=z1)}
	     z1})
	names(z)<-slnm
	assert.is(z,"list")
	class(z)<-c("list",paste(class(x),"list",sep="."))
	z
}

assert.is <- function(object, class2, text = ""){
	stopifnot(is.character(class2))
	if(missing(object))
		stop(paste(class2, "object missing in assert.is", text))
	stopifnot(length(class2) >= 1)

	if(!is_any(object = object, class2 = class2))
	{
		warning(paste("got ", class(object), sep = ""))
		message("got ", class(object), " when one of these classes was required:")
		print(class2)
		stop(text)
	}
}
assert.are <- function(object, class2, ...){
	assert.is(object, "list")
	if(!is_nothing(object)){
		for(obj in object)
			assert.is(object = obj, class2 = class2, ...)
	}
}
are <- function(object, class2){all(sapply(object, is_any, class2 = class2))}
are_null <- function(object){are(object, "NULL")}
is_vide<-is_nothing<-function(object){length(object)==0}
is_any <- function(object, class2){any(sapply(class2, is, object = object))}
setGeneric("is_prob", function(object,...) standardGeneric("is_prob"))
setMethod("is_prob",signature(object = "numeric" ), function(object,tolerance = 1e-3, ...){
         boo <- all(object >= 0 - tolerance & object <= 1 + tolerance, ...)
	 !is.na(boo) && boo})
are_prob <- function(object, ...){all(sapply(object, is_prob, ...))}
is_err <- function(object){is(object, "try-error")}
is_error <- function(object){is(object, "try-error")||is_nothing(object)}
is_unk <- function(object){all(is_any(object,c("try-error","NULL")))||all(is.na(object))||is_nothing(object)}
are_unk <- function(object){all(sapply(object, is_unk))}

nsize<-function(x){
	if(is_any(x,c("list","NULL","numeric","logical","character"))){c(1,length(x))}
	else if(is_any(x,c("matrix","array","data.frame"))){c(nrow(x),ncol(x))}
	else if(is_nothing(x)){c(0,0)}
	else {c(length(x),1)}
}

vect2string<-function(x,sep="",...){paste(x,sep="",collapse=sep,...)}
make_labels<-function(n,nmvar=c("X"),n.ini=1){n.ini<-n.ini[1]
    if(is_vide(n)){return(NULL)}
    stopifnot(n>0)
    paste(nmvar[1],c(1:n)+n.ini-1,sep="")
}

MakeNames<-function(x,nmvar=c("X","I"),...){
	k.MakeNames<-function(x,col=T,nmvar="I",force=FALSE,unique=T,n0=1){
		if(is_nothing(x)){return(x)}
		k.nms<-function(nm,nn){
			if(is_nothing(nm) | all(is_unk(nm)) | force==T ){
				nm<-make_labels(n=nn,nmvar=nmvar[1],n.ini=n0[1])}
			if(any(is.na(nm))){ind<-which(is.na(nm));nm[ind]<-paste(nmvar[1],ind,sep="")}
			if(unique==T){nm<-make.unique(nm)}
			nm	
		}
		if(is_any(x,c("numeric","logical","character","list"))){
			nm<-k.nms(nm=names(x),nn=length(x))
			names(x)<-nm;return(x)}
	
		assert.is(x,c("matrix","array","data.frame"))
		nn<-nsize(x)
		
		if(col){nm<-colnames(x);nn<-nn[2]}
		else{nm<-rownames(x);nn<-nn[1]}
	
		nm<-k.nms(nm=nm,nn=nn)
		
		if(col){colnames(x)<-nm}
		else{rownames(x)<-nm}
		x
}
	
	#-------------------------------------
	if(is_any(x,c("numeric","logical","character","list"))){
		x<-k.MakeNames(x=x,nmvar=nmvar[1],...)}
	else if(is_any(x,c("array","matrix","data.frame"))){
		x<-k.MakeNames(x,col=T,nmvar=nmvar[pmax(2,length(nmvar))],...)
		x<-k.MakeNames(x,col=F,nmvar=nmvar[1],...)
	}
	x
}
sameAsX_names<-function(x=NULL,y=NULL){
	assert.are(list(x,y),c("numeric","logical","character"))
	stopifnot(class(x)==class(y))
	
	x<-MakeNames(x,nmvar="X")
	y<-MakeNames(y,nmvar="X")
	stopifnot(length(x)>=length(y))
	if(are(list(x,y),"numeric")){zz<-as.numeric(NA)}
	else if(are(list(x,y),"character")){zz<-as.character(NA)}
	else if(are(list(x,y),"logical")){zz<-NA}
		
	z<-rep(zz,length(x));names(z)<-names(x)
	z[names(y)]<-y
	stopifnot(identical(names(x),names(z)))
	z
}

#---------------------------------
dabsTd<-function(x,df,ncp=0,...){dt(x=x,df=df,ncp=ncp,...)+dt(x=-x,df=df,ncp=ncp,...)}
attr(dabsTd,'name')<-'dabsTd'
setClass("est.lfdrmle", representation(LFDR.hat ="numeric",p0.hat = "numeric", ncp.hat= "numeric",
  stat="numeric",info="list"))
setValidity("est.lfdrmle", function(object) {
  dx<-object@LFDR.hat
  p0<-object@p0.hat
  ncp<-object@ncp.hat
  stx<-object@stat
  ix<-object@info
  nvar<-length(stx)
stopifnot(length(names(dx))>0)

stopifnot(all(c(length(p0),length(ncp))%in%c(1,nvar)))
stopifnot(length(stx)==length(dx)&&identical(names(stx),names(dx)))

if(!are_prob(list(dx[is.finite(dx)],p0[is.finite(p0)]))){stop("error:LFDR.hat or p0.hat out of [0,1]")}

})
new_est.lfdrmle <- function(LFDR.hat,p0.hat,ncp.hat,stat,method=NULL,info=list()){
	if(is_nothing(method)){method<-"lfdr.mle"}
	if(!is(info,"list")){info<-list()}
	assert.is(LFDR.hat,c("NULL","numeric","logical","try-error"))
	assert.is(stat,c("numeric","logical"))
	stopifnot(length(stat)>=length(LFDR.hat))
	
	if(is_error(LFDR.hat)){
		message(method," estimation failed ",Sys.time())
		LFDR.hat<-rep(as.numeric(NA),length(stat))
		names(LFDR.hat)<-names(stat)
		method<-paste("FAILED-",method,sep="")
		p0.hat<-ncp.hat<-as.numeric(NA)
		}
	
	nvar<-length(stat)


	if(is(stat,"logical")){xx<-as.numeric(stat);names(xx)<-names(stat);stat<-xx}
	if(is(LFDR.hat,"logical")){xx<-as.numeric(LFDR.hat);names(xx)<-names(LFDR.hat);LFDR.hat<-xx}
	if(is(p0.hat,"logical")){xx<-as.numeric(p0.hat);names(xx)<-names(p0.hat);p0.hat<-xx}
	if(is(ncp.hat,"logical")){xx<-as.numeric(ncp.hat);names(xx)<-names(ncp.hat);ncp.hat<-xx}	
	assert.are(list(LFDR.hat,p0.hat,ncp.hat,stat),"numeric")
	
  
	stat<-MakeNames(stat,nmvar="X")
	LFDR.hat<-MakeNames(LFDR.hat,nmvar="X")

	if(length(ncp.hat)==nvar){names(ncp.hat)<-names(LFDR.hat)}
	else if(length(ncp.hat)==1){names(ncp.hat)<-NULL}
	if(length(p0.hat)==nvar){names(p0.hat)<-names(LFDR.hat)}
	else if(length(p0.hat)==1){names(p0.hat)<-NULL}

	LFDR.hat<-sameAsX_names(y=LFDR.hat,x=stat)

	if(length(ncp.hat)>1&&length(stat)>1){
		ncp.hat<-sameAsX_names(y=ncp.hat,x=stat)
		p0.hat<-sameAsX_names(y=p0.hat,x=stat)
	stopifnot(identical(names(stat),names(p0.hat)) && identical(names(LFDR.hat),names(ncp.hat)))}
	stopifnot(all(c(length(p0.hat),length(ncp.hat))%in%c(1,nvar)))
	infox<-c(list(method=toupper(method)),info)
new("est.lfdrmle",LFDR.hat=LFDR.hat, p0.hat=p0.hat,ncp.hat=ncp.hat,stat=stat,info=infox)


}







###-----
setMethod("[", signature(x = "est.lfdrmle", i = "ANY", j = "missing"), function(x, i, j, drop){
	z<-x
	nvar<-length(x@LFDR.hat)
	z@LFDR.hat <- x@LFDR.hat[i];names(z@LFDR.hat)<-names(x@LFDR.hat)[i]
	z@stat <- x@stat [i];names(z@stat)<-names(x@LFDR.hat)[i]
	if(length(x@p0.hat)==nvar){z@p0.hat <- x@p0.hat[i];names(z@p0.hat)<-names(x@LFDR.hat)[i]}
	if(length(x@ncp.hat)==nvar){z@ncp.hat <- x@ncp.hat[i];names(z@ncp.hat)<-names(x@LFDR.hat)[i]}
	stopifnot(validObject(z))
	z})

setReplaceMethod("names", signature(x="est.lfdrmle",value="character"),function(x, value){
  nvar<-length(x@LFDR.hat)
  stopifnot(length(value)==nvar)
  names(x@LFDR.hat)<-value
  names(x@stat) <- value
  if(length(x@p0.hat)==nvar){names(x@p0.hat) <- value}
  if(length(x@ncp.hat)==nvar){names(x@ncp.hat) <- value}
  stopifnot(validObject(x))
	x})
                 
setMethod("names", signature(x = "est.lfdrmle"), function(x){names(x@LFDR.hat)})
setMethod("as.numeric", signature(x="est.lfdrmle"),function(x) {x@LFDR.hat})
setMethod("length", signature(x = "est.lfdrmle"), function(x){length(as.numeric(x))})

setMethod("is_prob",signature(object = "est.lfdrmle" ), function(object,tolerance = 1e-3, ...){
         is_prob(object@LFDR.hat,tolerance =tolerance,...) && is_prob(object@p0.hat,tolerance =tolerance,...)})

#based on Ye_new.r-------------------------------

k.log_wlik_mixture <- function(W,p0,dalt,dFUN,d0=0,w=1,...){sum(w*log(p0*dFUN(W,ncp=d0,...) +(1-p0)*dFUN(W,ncp=dalt,...)))}
k.get_wd_max <- function(W,p0,dFUN,lower.ncp=1/1e3, upper.ncp=100,d0=0,w=1,...){### get the maximum value of d_alt
 f <- function(dalt){k.log_wlik_mixture(p0=p0,W=W,dalt=dalt,w=w,dFUN=dFUN,d0=d0,...)}
 as.numeric(optimize(f, lower=lower.ncp, upper=upper.ncp,maximum=TRUE))[1]
}

k.get_wp0_max <- function(W,dalt,dFUN=dabsTd, lower.p0=0, upper.p0=1,d0=0,w=1, ...){### get the maximum value of p0
 f <- function(p0){k.log_wlik_mixture(p0=p0,W=W,dalt=dalt,w=w,dFUN=dFUN,d0=d0,...)}
 as.numeric(optimize(f, lower=lower.p0, upper=upper.p0,maximum=TRUE))[1]
}

k.get_wlfdr <- function(W,p0,dalt,dFUN,d0=0,w=1,...){
 zo<-(p0*dFUN(W,ncp=d0,...))/(p0*dFUN(W,ncp=d0,...)+(1-p0)*dFUN(W,ncp=dalt,...))^w
 names(zo)<-names(W);zo
}


k.N_wlik <- function( W, dFUN, lower.ncp, upper.ncp, lower.p0, upper.p0,d0,w=1,...){ 

	p0_max<-d_max<-NULL

	if(lower.p0==upper.p0){p0_max <- lower.p0}
	if(lower.ncp==upper.ncp){d_max <- lower.ncp}
	if(!is.null(p0_max)&&!is.null(d_max)){zo<-c(p0_max,d_max);names(zo)<-c("p0.max","ncp.max");return(zo)}

	if(is.null(p0_max)&&is.null(d_max)){		
		ff <- function(p0){
			dalt<-k.get_wd_max(p0=p0,W=W,w=w,dFUN=dFUN,lower.ncp=lower.ncp,upper.ncp=upper.ncp,d0=d0,...)
			k.log_wlik_mixture(p0=p0,W=W,w=w,dalt=dalt,dFUN=dFUN,d0=d0,...)}	
	p0_opt_result <- optimize(ff, lower=lower.p0, upper=upper.p0, maximum=TRUE)### p0 will  be in [lower.p0,upper.p0]
	p0_max <- as.numeric(p0_opt_result[1])
	}
	if(is.null(p0_max)&&!is.null(d_max)){		
		p0_max <- k.get_wp0_max(dalt=d_max, W=W,w=w,dFUN=dFUN, lower.p0=lower.p0, upper.p0=upper.p0,d0=d0,...)
	}
	if(!is.null(p0_max)&&is.null(d_max)){
		d_max <- k.get_wd_max(p0=p0_max, W=W,w=w,dFUN=dFUN, lower.ncp=lower.ncp, upper.ncp=upper.ncp,d0=d0,...)}
	### max.logL <- log_lik(p0_max, W, d_max,dFUN,df)
	if(!is_prob(p0_max)){stop("p0_max should be in [0,1]")}
	zo<-c(p0_max,d_max);#names(zo)<-c("p0.max","ncp.max");
	return(zo)
}





k.lfdr.mle<- function(stat,dFUN,lower.ncp=1/1e3,upper.ncp=200,lower.p0=0,upper.p0=1, d0=0,w=1,fixed.p0=NULL,fixed.ncp=NULL,...){#w: weights   
  
  
if(length(fixed.p0)>0){stopifnot(is_prob(fixed.p0));lower.p0<-upper.p0<-fixed.p0 }
if(length(fixed.ncp)>0){lower.ncp<-upper.ncp<-fixed.ncp}

  opt <- try(k.N_wlik(W=stat,w=w,d0=d0,dFUN=dFUN,lower.ncp=lower.ncp,upper.ncp=upper.ncp,lower.p0=lower.p0,upper.p0=upper.p0,...),silent=F)
  if(is_err(opt)){return(list())}
  lfdr_mix <- k.get_wlfdr(W=stat,w=1,p0=opt[1],dalt=opt[2],dFUN=dFUN,d0=d0,...)
	if(is_err(lfdr_mix)){return(list())}
  list(lfdr=lfdr_mix,p0=opt[1],ncp=opt[2])
}





#

k.perfeat<-function(stat,w=1,v,dFUN=dabsTd,d0,fixed.p0=NULL,fixed.ncp=NULL,lower.ncp=1/1e3, upper.ncp=200, lower.p0=0, upper.p0=1,...){
  
  assert.are(list(stat,d0,v,w), "numeric");nvars<-length(stat)
  stopifnot(is_prob(v))
  if(length(w)==1){w<-rep(w,length(stat))}
  if(abs(sum(w)-1)>1e-3){w<-w/sum(w)}
  wz<-wx<-w
  k.fun<-function(i){
    wx[i]<-v*wz[i];
    vars<- which(c(1:nvars)!=i)
    xt<-stat[vars];xp<-stat[i]

    stopifnot(all(wx>=0))
    zvar<-k.lfdr.mle(stat=stat,dFUN=dFUN,lower.ncp=lower.ncp,upper.ncp=upper.ncp,lower.p0=lower.p0,upper.p0=upper.p0,
           fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,d0=d0,w=wx,...)
    z1<-k.lfdr.mle(stat=xp,dFUN=dFUN,lower.ncp=lower.ncp,upper.ncp=upper.ncp,lower.p0=lower.p0,upper.p0=upper.p0,
           fixed.ncp=zvar$ncp,fixed.p0=zvar$p0,d0=d0,w=1,...)
    
    c(z1$lfdr,z1$p0,z1$ncp)
    
  }

  zo<-try(vapply(1:nvars,FUN=k.fun,FUN.VALUE=numeric(length=3)),silent=F)
  if(is_err(zo)){return(list())}
  colnames(zo)<-names(stat)
  
list(lfdr=zo[1,],p0=zo[2,],ncp=zo[3,])

}
k.lo<-function(stat,w=1,v,dFUN=dabsTd,d0,fixed.p0=NULL,fixed.ncp=NULL,lower.ncp=1/1e3, upper.ncp=200, lower.p0=0, upper.p0=1,...){
  z1<-k.lfdr.mle(stat = stat,d0=d0,w=w,dFUN=dFUN,fixed.p0=fixed.p0,lower.ncp=lower.ncp,
                upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,fixed.ncp=fixed.ncp,...)
  fixed.p0<-z1$p0
  assert.is(fixed.p0,"numeric");stopifnot(is_prob(fixed.p0))

  k.perfeat(stat=stat,w=w,v=v,dFUN=dFUN,d0=d0,fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,
      lower.ncp=lower.ncp,upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,...)
 
}

##--------
k.checkings<-function(stat,dFUN,lower.ncp,upper.ncp,lower.p0,upper.p0,fixed.p0=NULL,fixed.ncp=NULL,v=1,d0=0,w=1,...){
  if(all(is_unk(stat))){stop("invalid input")}
  assert.is(dFUN,"function")
 
  assert.are(list(stat,w,d0,v,lower.ncp,upper.ncp,lower.p0,upper.p0), "numeric")
  assert.are(list(fixed.p0,fixed.ncp), c("numeric","NULL"))

  nvars<-length(stat);if(nvars==0){return(NULL)}
  stat<-MakeNames(stat,nmvar="X")
  
  if(length(w)==1){w<-rep(w,nvars)}
  stopifnot(length(w)%in%c(1,nvars))
  stopifnot(all(w>=0))

  W<-nW<-stat

  pW<-dFUN(stat,ncp=0,...)
  ind<-which(is.finite(stat) & is.finite(pW))
  if(length(ind)==0){stop("non finite statistic values")}
  stat<-stat[ind];w<-w[ind]

  if(length(fixed.p0)>0){stopifnot(is_prob(fixed.p0));lower.p0<-upper.p0<-fixed.p0 }
  if(length(fixed.ncp)>0){lower.ncp<-upper.ncp<-fixed.ncp}
  if(!are_prob(list(lower.p0, upper.p0,v))){stop("v, lower.p0 and upper.p0 should be in [0,1]")}
  stopifnot(lower.p0 <= upper.p0 && lower.ncp <= upper.ncp)
 
  nvars<-length(stat)
  stopifnot(length(w)%in%c(nvars))
  stopifnot(all(w>=0))
  if(abs(sum(w)-1)>1e-3){w<-w/sum(w)}
  if(abs(sum(w)-1)>1e-3){stop("problem with weights")}
    
  siz.num1<-(sapply(list(d0,v,lower.ncp,upper.ncp,lower.p0,upper.p0),length))
  names(siz.num1)<-c("d0","v","lower.ncp","upper.ncp","lower.p0","upper.p0")
  siz.num0<-(sapply(list(fixed.p0,fixed.ncp),length))
  names(siz.num0)<-c("fixed.p0","fixed.ncp")
  siz.numn<-(sapply(list(stat,w),length))
  names(siz.numn)<-c("stat","w")
  if(!all(siz.num1%in%c(1))){
    stop(paste("numeric inputs: ",vect2string(names(siz.num1[!siz.num1%in%1]),sep=" ")," have length !=1"))}
  if(!all(siz.num0%in%c(0,1))){
    stop(paste("numeric inputs: ",vect2string(names(siz.num0[!siz.num0%in%c(0,1)]),sep=" ")," have length !=0 OR 1"))}
  if(!all(siz.numn%in%c(nvars))){
    stop(paste("numeric inputs: ",vect2string(names(siz.num0[!siz.num0%in%c(nvars)]),sep=" ")," have length !=",nvars))}
         
    list(stat=stat,orig.stat=W,dFUN=dFUN,w=w,lower.ncp=lower.ncp,upper.ncp=upper.ncp,lower.p0=lower.p0,upper.p0=upper.p0,
         fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,v=v,d0=d0)
}


##---------------
lfdr.hats<-function(stat=NULL,lfdr.fun="lfdr.mle",dFUN=dabsTd,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,
		    fixed.p0=NULL,fixed.ncp=NULL,d0=0,v=1,w=1,...){
	
	z<-k.checkings(stat=stat,w=w,v=v,dFUN=dFUN,d0=d0,lower.ncp=lower.ncp,upper.ncp=upper.ncp,
                 lower.p0=lower.p0, upper.p0= upper.p0,fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,...)

  #---------------
	if(lfdr.fun%in%c("lfdr.mle","mle","L0O","l0o")){z$w<-1
	zo<-k.lfdr.mle(stat=z$stat,dFUN=z$dFUN,lower.ncp=z$lower.ncp,upper.ncp=z$upper.ncp,lower.p0=z$lower.p0,
		   upper.p0=z$upper.p0,fixed.p0=z$fixed.p0,fixed.ncp=z$fixed.ncp,d0=z$d0,w=z$w,...)
	method<-"MLE"
	}
	else if(lfdr.fun%in%c("lfdr.mdl","mdl")){z$w<-1;z$v<-0
	  zo<-k.perfeat(stat=z$stat,w=z$w,v=z$v,dFUN=z$dFUN,d0=z$d0,fixed.p0=z$fixed.p0,fixed.ncp=z$fixed.ncp,
	      lower.ncp=z$lower.ncp,upper.ncp=z$upper.ncp, lower.p0=z$lower.p0, upper.p0= z$upper.p0,...)
	       method<-"MDL"}
      
	else if(lfdr.fun%in%c("lfdr.l1o","l1o","L1O")){z$w<-1;z$v<-0
	zo<-k.lo(stat=z$stat,w=z$w,v=z$v,dFUN=z$dFUN,d0=z$d0,fixed.p0=z$fixed.p0,fixed.ncp=z$fixed.ncp,
	      lower.ncp=z$lower.ncp,upper.ncp=z$upper.ncp, lower.p0=z$lower.p0, upper.p0= z$upper.p0,...)
	method<-"L1O"
	}
	else if(lfdr.fun%in%c("lfdr.lho","lho","LHO")){z$w<-1;z$v<-1/2
	zo<-k.lo(stat=z$stat,w=z$w,v=z$v,dFUN=z$dFUN,d0=z$d0,fixed.p0=z$fixed.p0,fixed.ncp=z$fixed.ncp,
	      lower.ncp=z$lower.ncp,upper.ncp=z$upper.ncp, lower.p0=z$lower.p0, upper.p0= z$upper.p0,...)
	method<-"LHO"}
	else if(lfdr.fun%in%c("lfdr.lo","lo","LO")){#z$w<-1;z$v<-0
	zo<-k.lo(stat=z$stat,w=z$w,v=z$v,dFUN=z$dFUN,d0=z$d0,fixed.p0=z$fixed.p0,fixed.ncp=z$fixed.ncp,
	      lower.ncp=z$lower.ncp,upper.ncp=z$upper.ncp, lower.p0=z$lower.p0, upper.p0= z$upper.p0,...)
	method<-paste("LO",round(z$v,1),sep="-v")
	}
      else if(lfdr.fun%in%c("lfdr.mdlo","mdlo","MDLO")){#z$w<-1;z$v<-0
	zo<-k.perfeat(stat=z$stat,w=z$w,v=z$v,dFUN=z$dFUN,d0=z$d0,fixed.p0=z$fixed.p0,fixed.ncp=z$fixed.ncp,
	      lower.ncp=z$lower.ncp,upper.ncp=z$upper.ncp, lower.p0=z$lower.p0, upper.p0= z$upper.p0,...)
	method<-paste("MDLO",round(z$v,1),sep="-v")
	}	

	info<-z[c("lower.ncp","upper.ncp","lower.p0","upper.p0","fixed.p0","fixed.ncp","d0","v")]#,"w"
	lix<-sapply(info,length)
	info<-info[lix>0]

    
  #--------------------
 zo<-new_est.lfdrmle(LFDR.hat=zo$lfdr,p0.hat=zo$p0,ncp.hat=zo$ncp,stat=z$orig.stat,method=method,info=info)
 est2list(zo)
}




#---------------
lfdr.mle<-function(x,dFUN=dabsTd, lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.p0=NULL,fixed.ncp=NULL,d0=0,...){
  lfdr.hats(lfdr.fun="lfdr.mle",stat=x,dFUN=dFUN,lower.ncp=lower.ncp,upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,
            fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,d0=d0,w=1,v=1,...)
}
lfdr.mdl<-function(x,dFUN=dabsTd,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.p0=NULL,fixed.ncp=NULL,d0=0,...){
	lfdr.hats(lfdr.fun="lfdr.mdl",stat=x,dFUN=dFUN, w=1,v=0,d0=d0,lower.ncp=lower.ncp,upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,
		fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,...)
}
lfdr.l1o<-function(x,dFUN=dabsTd,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.p0=NULL,fixed.ncp=NULL,d0=0,...){
  lfdr.hats(lfdr.fun="lfdr.l1o",stat=x,dFUN=dFUN,w=1,v=0,d0=d0,lower.ncp=lower.ncp,upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,
            fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,...)
}
lfdr.lho<-function(x,dFUN=dabsTd,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.p0=NULL,fixed.ncp=NULL,d0=0,...){
  lfdr.hats(lfdr.fun="lfdr.lho",stat=x,dFUN=dFUN,w=1,v=1/2,d0=d0,lower.ncp=lower.ncp,upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,
            fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,...)
}
lfdr.lo<-function(x,dFUN=dabsTd,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.p0=NULL,fixed.ncp=NULL,v=0,d0=0,...){
  lfdr.hats(lfdr.fun="lfdr.lo",stat=x,dFUN=dFUN,w=1,v=v,d0=d0,lower.ncp=lower.ncp,upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,
            fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,...)
}
lfdr.mdlo<-function(x,v=0,dFUN=dabsTd,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.p0=NULL,fixed.ncp=NULL,d0=0,...){
  lfdr.hats(lfdr.fun="lfdr.mdlo",stat=x,dFUN=dFUN,w=1,v=v,d0=d0,lower.ncp=lower.ncp,upper.ncp=upper.ncp, lower.p0=lower.p0, upper.p0= upper.p0,
            fixed.p0=fixed.p0,fixed.ncp=fixed.ncp,...)
}

attr(lfdr.l1o,"name")<-"L1O"
attr(lfdr.mdl,"name")<-"MDL"
attr(lfdr.mle,"name")<-"MLE"
attr(lfdr.lho,"name")<-"LHO"
#--from Ye### general function 
vectorized.dabsTd <- function(ncp, ...) {sapply(ncp, function(ncp) {dabsTd(ncp = ncp, ...)})}
get_n.groups <- function(n.features,n.null,group.size){#### get group size 
	n <- (n.features - n.null)
	if(group.size >n){
	out <- 1
	}else{
	if((n %% group.size)!=0) {
	  out <- floor(n/group.size)+1
	 }else{
	  out <- n/group.size
	 }
	}
       out
}

get_groups <- function(n.features,n.null,n.groups){### evenly divided groups
	n <- (n.features - n.null)
	out <- rep(floor(n / n.groups), n.groups)
	if((n %% n.groups) != 0) out[1:(n %% n.groups)] <- out[1:(n %% n.groups)] + 1
	out
}

 
get_null.n.groups <- function(n.null,null_group.size){### get null groups
	n <- n.null
	if(null_group.size >n){
	out <-1
	}else{
	if((n %% null_group.size)!=0){
	 out <- floor(n/null_group.size)+1
	  }else{
	   out <- n/null_group.size
	  }
	 }
	out
}



get_null.groups <- function(n.null,null.n.groups){### get null groups list
	n <- n.null
	out <- rep(floor(n / null.n.groups), null.n.groups)
	if((n %% null.n.groups) != 0) out[1:(n %% null.n.groups)] <- out[1:(n %% null.n.groups)] + 1
	out
}


rFUN_generator <- function(base_rFUN,...){  
	function(n, ncp,...)
	{
   stopifnot(length(n) == length(ncp))
    unlist(lapply(1:length(n), function(i) base_rFUN(n[i],  ncp[i],...)))  
 }
}

simulated_W <- function(n.features,n.null,rFUN=rFUN_generator(rchisq),true.ncp1,df, #can deal with different true.ncp1 
	N=(n.features-n.null),sided=ifelse(identical(rFUN, rFUN_generator(rchisq)),1,2),logic=TRUE){ 
	  if(n.null != 0 & logic == TRUE)
	  {
	    true.ncp1 <- c(rep(0,(length(N)-length(true.ncp1))), true.ncp1)
	   #### N <- c(n.null, N)
	  }
	  if(logic)
	  {
	  W <- numeric(n.features)
	  temp <- rFUN(rep(1, length(N)), df=df, true.ncp1)
	  W <- unlist(lapply(1:length(N), function(i) rep(temp[i], N[i])))
	  } else {
	    ### W <- numeric(n.null + n.features)
	  W <- numeric(n.features)
	  if(n.null != 0) W[1:n.null] <- rFUN(n.null, df, 0)
	  temp <- rFUN(rep(1, length(N)), df, true.ncp1)
	  W[(n.null+1):n.features] <- unlist(lapply(1:length(N), function(i) rep(temp[i], N[i])))
	  }
	  if(sided==1)
	  {
	    Data <- W
	  } else {
	    Data <- abs(W)
	  }
	  return(Data)
}



get_simulated_W  <- function(n.features,true.p0,list.ncp1,rFUN,df){
#### simulate indicator
	II <- ifelse(runif(n.features,min=0,max=1) < true.p0, 0,1)
	
      #### simulate true.ncp1
	temp <- sample(list.ncp1,size =length(II),replace =TRUE)
	true.ncp1 <- ifelse(II==0,rep(0,length(II)), temp)
	
      #### simulated W
	W<- unlist(lapply(1:n.features, function(i) rFUN(n=1, df = df, ncp = true.ncp1[i])))
	return(cbind(W,II,true.ncp1))
	
}



#from Ye_new.r-------------------------
mix.ye.mle<-function(W,dFUN,df,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.ncp=NULL,fixed.p0=NULL){
	tol<-1e-3
	if(!is.null(fixed.p0)){lower.p0<-fixed.p0-tol/2;upper.p0<-fixed.p0+tol/2}
	if(!is.null(fixed.ncp)){lower.ncp<-fixed.ncp-tol/2;upper.ncp<-fixed.ncp+tol/2}
	lower<-lower.ncp;upper<-upper.ncp
	    #W<-MakeNames(W);W<-W[is.finite(W)]
	ostat<-W
	z<-k.checkings(stat=W,dFUN=dFUN,lower.ncp=lower,upper.ncp=upper,lower.p0=lower.p0,upper.p0=upper.p0,df=df)
	W<-z$stat
	#=============
	new.estimation <- function(LFDR.hat,p0.hat,ncp1.hat,d.hat){
		list(LFDR.hat=LFDR.hat, p0.hat=p0.hat,ncp.hat=ncp1.hat)
		#new("estimation",LFDR.hat=LFDR.hat, p0.hat=p0.hat,ncp1.hat=ncp1.hat,d.hat=d.hat)
	}
	log_lik_mixture <- function(p0,W,d_alt,dFUN,df=1){
		sum(log(p0*dFUN(W,df=df,ncp=0) + (1-p0)*dFUN(W,df=df,ncp=d_alt)))
	}
	  
	### get the maximum value of d_alt
	get_d_max <- function(p0,W,dFUN, df=1,lower,upper,...){
		f <- function(d_alt){log_lik_mixture(p0=p0,W=W,d_alt=d_alt,dFUN=dFUN,df=df)}
		as.numeric(optimize(f, lower=lower, upper=upper,maximum=TRUE,...))[1]
	}
	  
	  
	### get the maximum value of p0
	N_lik <- function(W,dFUN,df=1,lower,upper,lower.p0,upper.p0,...){
		ff <- function(p0){
			log_lik_mixture(p0,W,d_alt=get_d_max(p0,W=W,dFUN=dFUN,df=df,lower=lower,upper=upper,...),
					dFUN=dFUN,df=df)
		}
		### the value of p0 can only be in the interval [lower.p0,upper.p0]
		p0_opt_result <- optimize(ff, lower=lower.p0, upper=upper.p0, maximum=TRUE,...)
		p0_max <- as.numeric(p0_opt_result[1])
		d_max <- get_d_max(p0_max, W=W,dFUN=dFUN,df=df,lower=lower, upper=upper, ...)
		     ### max.logL <- log_lik(p0_max, W, d_max,dFUN,df)
		return(c(p0_max,d_max))
	}
	  
	###### calculate the local false discovery rate
	  
	get_lfdr <- function(pval,qFUN,W,p0,d_alt,dFUN=dchisq,df=1){
		if(missing(W) && is.function(qFUN)){
			W <- qFUN(pval, df = df,lower.tail=FALSE)
		}
		assert.is(W, "numeric")
	       (p0*dFUN(W,df=df,ncp=0))/(p0*dFUN(W,df=df,ncp=0)+(1-p0)*dFUN(W,df=df,ncp=d_alt))
	}
	  
	  
	##### the function of all the functions above
	MLE.LFDR_mixture <- function(pval,qFUN,W,df=1,dFUN=dchisq,lower,upper,ower.p0,upper.p0,...){
		if(missing(W) && is.function(qFUN))
			W <- qFUN(pval, df = df,lower.tail=FALSE)
		else
			stopifnot(missing(pval) && missing(qFUN))
		assert.is(W, "numeric")
		opt <- N_lik(W=W,dFUN=dFUN,df=df,lower=lower,upper=upper,lower.p0=lower.p0,upper.p0=upper.p0)
		lfdr_mix <- get_lfdr(W=W,p0=opt[1],d_alt=opt[2],dFUN=dFUN,df=df)
		return(list(lfdr_mix,opt[1],opt[2]))
	}
	  
	#### mixture model
	mixture.mle.estimation <- function(pval,qFUN,W,df,dFUN,lower,upper,lower.p0,upper.p0,push,...){                    
		if(missing(W) && is.function(qFUN)) {W <- qFUN(p=pval, df = df,ncp=0,lower.tail=FALSE)}
		assert.is(W, "numeric")
		opt <- N_lik(W=W,dFUN=dFUN,df=df,lower=lower,upper=upper,lower.p0=lower.p0,upper.p0=upper.p0,...)
		lfdr_mix <- get_lfdr(W=W,p0=opt[1],d_alt=opt[2],dFUN=dFUN,df=df)
		d <- numeric(0)
		out <- new.estimation(LFDR.hat=lfdr_mix,p0.hat=opt[1],ncp1.hat=opt[2],d.hat=d)
		if(!missing(push)) push(list(out=out, W=W))
		out  
	  }
	#------------------------
	Wx<-dFUN(W,df=df)
	W<-W[is.finite(Wx)]
	    
	opt <- N_lik(W=W,dFUN=dFUN,df=df,lower=lower,upper=upper,lower.p0=lower.p0,upper.p0=upper.p0)
	lfdr_mix <- get_lfdr(W=W,p0=opt[1],d_alt=opt[2],dFUN=dFUN,df=df)
	list(LFDR.hat=sameAsX_names(y=lfdr_mix,x=ostat),p0.hat=opt[1],ncp.hat=opt[2]) 
	    
}
pure.ye.mle<-function(W,dFUN,df,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,fixed.ncp=NULL,fixed.p0=NULL){

	tol<-1e-3
	if(!is.null(fixed.p0)){lower.p0<-fixed.p0-tol/2;upper.p0<-fixed.p0+tol/2}
	if(!is.null(fixed.ncp)){lower.ncp<-fixed.ncp-tol/2;upper.ncp<-fixed.ncp+tol/2}
	lower<-lower.ncp;upper<-upper.ncp
	  #W<-MakeNames(W);W<-W[is.finite(W)]
	ostat<-W
	z<-k.checkings(stat=W,dFUN=dFUN,lower.ncp=lower,upper.ncp=upper,lower.p0=lower.p0,upper.p0=upper.p0,df=df)
	W<-z$stat
	#=============
	log_lik_mixture <- function(p0,W,d_alt,dFUN,df=1){
		sum(log(p0*dFUN(W,df=df,ncp=0) + (1-p0)*dFUN(W,df=df,ncp=d_alt)))
	}
      
	### get the maximum value of d_alt
	get_d_max <- function(p0,W,dFUN, df=1,lower,upper,...){
		f <- function(d_alt){log_lik_mixture(p0=p0,W=W,d_alt=d_alt,dFUN=dFUN,df=df)}
		as.numeric(optimize(f, lower=lower, upper=upper,maximum=TRUE,...))[1]
	}
            
	### get the maximum value of p0
	N_lik <- function(W,dFUN,df=1,lower,upper,lower.p0,upper.p0,...){
		ff <- function(p0){
			log_lik_mixture(p0,W,d_alt=get_d_max(p0,W=W,dFUN=dFUN,df=df,lower=lower,upper=upper,...),
					dFUN=dFUN,df=df)
			}
		### the value of p0 can only be in the interval [lower.p0,upper.p0]
		p0_opt_result <- optimize(ff, lower=lower.p0, upper=upper.p0, maximum=TRUE,...)
		p0_max <- as.numeric(p0_opt_result[1])
		d_max <- get_d_max(p0_max, W=W,dFUN=dFUN,df=df,lower=lower, upper=upper, ...)
		### max.logL <- log_lik(p0_max, W, d_max,dFUN,df)
		return(c(p0_max,d_max))
	}
	###### calculate the local false discovery rate
	
	get_lfdr <- function(pval,qFUN,W,p0,d_alt,dFUN=dchisq,df=1){
		if(missing(W) && is.function(qFUN)){
			W <- qFUN(pval, df = df,lower.tail=FALSE)
		}
		assert.is(W, "numeric")
		(p0*dFUN(W,df=df,ncp=0))/(p0*dFUN(W,df=df,ncp=0)+(1-p0)*dFUN(W,df=df,ncp=d_alt))
	}
           
	##### the function of all the functions above
	MLE.LFDR_mixture <- function(pval,qFUN,W,df=1,dFUN=dchisq,lower,upper,lower.p0,upper.p0,...){
		if(missing(W) && is.function(qFUN))
			W <- qFUN(pval, df = df,lower.tail=FALSE)
		else
			stopifnot(missing(pval) && missing(qFUN))
		assert.is(W, "numeric")
		opt <- N_lik(W=W,dFUN=dFUN,df=df,lower=lower,upper=upper,lower.p0=lower.p0,upper.p0=upper.p0)
		lfdr_mix <- get_lfdr(W=W,p0=opt[1],d_alt=opt[2],dFUN=dFUN,df=df)
		return(list(lfdr_mix,opt[1],opt[2]))
	}
	#-----------------------
	Wx<-dFUN(W,df=df)
	W<-W[is.finite(Wx)]
	opt <- N_lik(W=W,dFUN=dFUN,df=df,lower=lower,upper=upper,lower.p0=lower.p0,upper.p0=upper.p0)
	lfdr_mix <- get_lfdr(W=W,p0=opt[1],d_alt=opt[2],dFUN=dFUN,df=df)
	
	return(list(LFDR.hat=sameAsX_names(x=ostat,y=lfdr_mix),p0.hat=opt[1],ncp.hat=opt[2]))
}
ye.mle<-function(W,df,dFUN=dabsTd,lower.ncp=1e-3,upper.ncp=20,lower.p0=0,upper.p0=1,d0=0,fixed.ncp=NULL,fixed.p0=NULL,...){
	tol<-1e-3
	if(!is.null(fixed.p0)){lower.p0<-fixed.p0-tol/2;upper.p0<-fixed.p0+tol/2}
	if(!is.null(fixed.ncp)){lower.ncp<-fixed.ncp-tol/2;upper.ncp<-fixed.ncp+tol/2}
	#W<-MakeNames(W);W<-W[is.finite(W)]
	ostat<-W
	z<-k.checkings(stat=W,dFUN=dFUN,lower.ncp=lower.ncp,upper.ncp=upper.ncp,lower.p0=lower.p0,upper.p0=upper.p0,df=df)
	W<-z$stat
	#=============
	log_lik_mixture <- function(p0,W,dalt,dFUN,df,d0=0,...){
		sum(log(p0*dFUN(W,df=df,ncp=d0,...) + (1-p0)*dFUN(W,df=df,ncp=dalt,...)))
	}
	
	get_d_max <- function(p0,W,dFUN, df,lower.ncp,upper.ncp,d0=0,...){### get the maximum value of d_alt
		f <- function(dalt){log_lik_mixture(p0=p0,W=W,dalt=dalt,dFUN=dFUN,df=df,d0=d0,...)}
		as.numeric(optimize(f, lower=lower.ncp, upper=upper.ncp,maximum=TRUE))[1]
	}
	
	N_lik <- function(W,dFUN,df,lower.ncp,upper.ncp,lower.p0,upper.p0,d0=0,...){### get the maximum value of p0
		ff <- function(p0){
			log_lik_mixture(p0,W,dalt=get_d_max(p0,W=W,dFUN=dFUN,df=df,lower.ncp=lower.ncp,upper.ncp=upper.ncp,d0=d0,...),
					dFUN=dFUN,df=df,d0=d0,...)}
		### the value of p0 can only be in the interval [lower.p0,upper.p0]
		p0_opt_result <- optimize(ff, lower=lower.p0, upper=upper.p0, maximum=TRUE)
		p0_max <- as.numeric(p0_opt_result[1])
		d_max <- get_d_max(p0=p0_max, W=W,dFUN=dFUN,df=df,lower.ncp=lower.ncp, upper.ncp=upper.ncp,d0=d0, ...)
		### max.logL <- log_lik(p0_max, W, d_max,dFUN,df)
		return(c(p0_max,d_max))
	}
	
	get_lfdr <- function(pval,qFUN,W,p0,dalt,dFUN,df,d0,...){###### calculate the local false discovery rate
		if(missing(W) && is.function(qFUN)){W <- qFUN(pval=pval, df = df,lower.tail=FALSE,...)}
		assert.is(W, "numeric")
		(p0*dFUN(W,df=df,ncp=d0,...))/(p0*dFUN(W,df=df,ncp=d0,...)+(1-p0)* dFUN(W,df=df,ncp=dalt,...))
	}
	
	
	mixture.mle.estimation <- function(W,df,dFUN,lower.ncp,upper.ncp,lower.p0,upper.p0,d0,...){                    
		assert.is(W, "numeric")
		opt <- N_lik(W=W,dFUN=dFUN,df=df,lower.ncp=lower.ncp,upper.ncp=upper.ncp,lower.p0=lower.p0,upper.p0=upper.p0,d0=d0,...)
		lfdr_mix <- get_lfdr(W=W,p0=opt[1],dalt=opt[2],dFUN=dFUN,df=df,d0=d0,...)
		d <- numeric(0)    
		infox<-list(method="ye.mle")
		new_est.lfdrmle(stat=ostat,LFDR.hat=lfdr_mix,p0.hat=opt[1],ncp.hat=opt[2],info=infox)
	}
  
	zo<-mixture.mle.estimation (W=W,df=df,dFUN=dFUN,lower.ncp=lower.ncp,upper.ncp=upper.ncp,lower.p0=lower.p0,upper.p0=upper.p0,d0=d0,...)
	est2list(zo)}

#----------

