 print.info<-echo<-FALSE

#-------------classes 1d,2d------
setClassUnion("d1classx",c("numeric","character","logical","integer",'Numeric'))
setClassUnion("d2classx",c("matrix","data.frame",'Matrix'))
d1class<-c("numeric","character","logical","integer",'Numeric')
d2class<-c("matrix","data.frame",'Matrix')
is.d1class<-function(x){Is(x,d1class)}
is.d2class<-function(x){Is(x,d2class)}
#-------------------
setMethod('nrow', signature(x = "d1classx" ), function(x){1})#length(x)
setMethod('ncol', signature(x = "d1classx" ), function(x){length(x)})
setMethod('nrow', signature(x = "NULL" ), function(x){0})
setMethod('ncol', signature(x = "NULL" ), function(x){0})
setMethod("rownames", signature(x = "d1classx"), function(x){NULL})#names(x)
setMethod("colnames", signature(x = "d1classx"), function(x){names(x)})

setReplaceMethod("colnames", signature(x="d1classx",value="d1classx"),function(x, value) {names(x)<-value})
#

setGeneric("fnames", function(object, value) standardGeneric("fnames"))
setGeneric("fnames<-", function(object, value) standardGeneric("fnames<-"))		
setGeneric("pnames", function(object, value) standardGeneric("pnames"))
setGeneric("pnames<-", function(object, value) standardGeneric("pnames<-"))
setMethod("pnames", signature(object = "d2classx"), function(object){colnames(object)})
setMethod("fnames", signature(object = "d2classx"), function(object){rownames(object)})
setMethod("pnames", signature(object = "d1classx"), function(object){names(object)})
setMethod("fnames", signature(object = "d1classx"), function(object){NULL})
setReplaceMethod("fnames", signature(object="d2classx"),function(object, value) {rownames(object)<-value})
setReplaceMethod("pnames", signature(object="d2classx"),function(object, value) {colnames(object)<-value})

#nvar<-function(x){}
setGeneric("nvar", function(x) standardGeneric("nvar"))
setMethod("nvar", signature(x="d1classx"),function(x) {1})
setMethod("nvar", signature(x="d2classx"),function(x) {nrow(x)})

#nsam<-function(x){}
setGeneric("nsam", function(x) standardGeneric("nsam"))
setMethod("nsam", signature(x="d1classx"),function(x) {length(x)})
setMethod("nsam", signature(x="d2classx"),function(x) {ncol(x)})
#is.paired<-function(x,y){}
#is.positive<-function(x){}
setGeneric("is.paired", function(x,y) standardGeneric("is.paired"))
setGeneric("is.positive", function(x) standardGeneric("is.positive"))
setMethod("is.positive", signature(x = "d2classx"), function(x){
	x<-as.numeric(x)
	is.positive(x)})
setMethod("is.positive", signature(x = "d1classx"), function(x){
	all(x[is.finite(x)]>0)})


#

nsize<-sizex<-function(x){
	if(is.d1class(x) | Is(x,c('list','NULL'))){c(1,length(x))}
	else if(is.d2class(x)){c(nrow(x),ncol(x))}
	else if(is.vide(x)){c(0,0)}
	else {c(length(x),1)}
}
nelem<-function(x,...){nn<-nsize(x);nn[1]*nn[2]}
between<-entre<-function(x,lower,upper){
    assert.are(list(lower,upper),'numeric')
    assert.is(x,c('numeric','logical'))
    #if(!all(is.finite(x))){return(NA)}
    z<-sort(c(lower,upper))
    
    x >= lower & x <= upper
}




#-------------------------------------
grep.or<-function(x,pattern,fixed=FALSE,exact=FALSE,ind=T,unik=T,...){
	z<-lapply(1:length(pattern),FUN=function(i){
	  if(exact){indx<-which(x%in%pattern[i])}
	  else{indx<-grep(x=x,pattern=pattern[i],fixed=fixed,...)}
	  if(!ind){x[indx]} else {indx}})
	zo<-unlist(z)
	if(unik){zo<-unique(zo)}
		zo
}
grepl.or<-function(x,pattern,fixed=FALSE,exact=FALSE,unik=T,...){z<-grep.or(x=x,pattern=pattern,fixed=fixed,exact=exact,ind=T,unik=T,...);length(z)>0}
string2vect<-function (x, sep = ",", fixed  = TRUE) {	
    assert.is(x, d1class)
    if(!is(x,'character')){x<-as(x,'chararcter')}
	zo<-strsplit(x=x, split=sep, fixed = fixed , perl = FALSE, useBytes = FALSE)
    sapply(1:length(zo),FUN=function(i){zo[[i]]})
}
vect2string<-function(x,sep=''){
	nx<-length(x);zo<-''
	for (i in 1:nx){
		if(i==1){zo<-x[i]}
		else{zo<-paste(zo,x[i],sep=sep)}
	}
zo}
#CorrectLimits<-function(x){}
setGeneric("CorrectLimits", function(x,...) standardGeneric("CorrectLimits"))
setMethod("CorrectLimits", signature(x = "numeric" ), function(x,maxlim=1,minlim=0){
	 stopifnot(!are.null(list(maxlim,minlim)))
	x[x<minlim]<-minlim
	x[x>maxlim]<-maxlim
	x})


k.rm.elem<-function(x,nm.rem='',row=T){
    if(is.vide(nm.rem)){return(x)}
    if(is(nm.rem,'numeric')){indx<-nm.rem}
    if(Is(x,c('list',d1class))){
        if(is(nm.rem,'character')){indx<-which(names(x)%in%nm.rem)}
        if(length(indx)>0){x<-x[-indx]}
        return(x)    
    }
    else if(Is(x,c(d2class)) && row){
        if(is(nm.rem,'character')){indx<-which(rownames(x)%in%nm.rem)}
        if(length(indx)>0){x<-x[-indx,,drop=FALSE]}
        return(x)   
    }
    else if(Is(x,c(d2class)) && !row){
        if(is(nm.rem,'character')){indx<-which(colnames(x)%in%nm.rem)}
        if(length(indx)>0){x<-x[,-indx,drop=FALSE]}
        return(x)   
    }
}
k.uniquex<-function(x,y=NULL,vip=1){assert.is(vip,'numeric')
	nm<-names(x);unm<-unique(nm)
	if(is.unique(names(x))){return(x)}
	lnm<-table(nm)
	z<-lapply(1:length(unm),FUN=function(i){
		indx<-which(nm%in%unm[i]);ix<-min(c(vip,length(indx)))
		#if(length(indx)>1){ix<-min(c(vip,length(indx)))}
		#else if(length(indx)==1){ix<-1}
		x[[indx[ix]]]})
	names(z)<-unm
         z}

nunique<-uniquex<-function(x,y=NULL,vip=1){x<-c(x,y)
	if(Is(x,d1class)){zo<-unique(x)}
	else if(is.vide(x)){zo<-x}	
	else if(is(x,'list')){zo<-k.uniquex(x=x)}
	assert.are(list(x),c('list','data.frame','NULL'))
	zo
}






#
get.elem<-function(x,nm='info'){nm<-nm[1]
  assert.is(x,'list');assert.is(nm,'character');zo<-NULL
  if (is(x,'list')&&nm%in%names(x)){eval(parse(text=paste('zo<-x$',nm,sep='')))}
  else if (isS4(x)){  #!Is(x,c(d1class,d2class,'list'))&&
	if(nm%in%slotNames(x)){eval(parse(text=paste('zo<-x@',nm,sep='')))}
        else{stop(paste('get.elem: bad nm (',nm,') or x'),sep='')}}
  else{stop('get.elem: bad x')}#XXX|:no visible binding for global variable ÔzoÕ
  zo
}
get.info<-function(x){get.elem(x=x,nm='info')}
get.lfdr<-function(x){get.elem(x=x,nm='LFDR.hat')}
get.p0<-function(x){get.elem(x=x,nm='p0.hat')}
get.ncp<-function(x){get.elem(x=x,nm='ncp.hat')}
#get.param<-function(x){}
setGeneric("get.param", function(x) standardGeneric("get.param"))
sameAsY<-function(x,y){
	assert.are(list(x,y),c(d1class,'NULL'))
	if(length(y)==0){stop('sameAsY: empty y')}
	if(length(x)==0){stop('sameAsY: empty x')}
	x<-MakeNames(x);y<-MakeNames(y)
	#z<-k.indSort(x=x,ref=y,ind=FALSE, inter=FALSE)
	z<-k.indSort(x=names(x),ref=names(y),ind=T, inter=FALSE)

	xx<-x[z$x];yy<-y[z$ref]
	names(xx)<-names(yy)<-z$names

	ok<-identical(names(xx),names(yy))
	stopifnot(ok)
	xx
}


##----Dvid:------
default<-function (object, name, verbose, return.value = object) {
    if (missing(verbose)) 
        verbose <- FALSE
    if (verbose) {
        if (missing(name)) 
            name <- "actual argument"
        if (is.function(name)) 
            name <- name(object)
        stopifnot(length(name) == 1 && is.character(name))
        prefix <- paste(name, "was set to default value of ")
        suffix <- paste(" on ", date(), ".", sep = "")
        object.name <- if (isS4(object) && any(c("annotation", 
            "ann") %in% slotNames(object))) {
            object.nam <- try(paste("object of annotation '", 
                annotation(object), "'", sep = ""))
            if (!is(object.nam, "character")) {
                message("bad object.nam from object0")
                browser()
            }
            object.nam
        }
        else object
        if (length(object.name) == 1 && (is(object.name, "character") || 
            is(object.name, "numeric"))) 
            message(prefix, object.name, suffix)
        else {
            cat("\n", prefix)
            print(object.name)
            cat(paste(suffix, "\n\n"))
        }
    }
    return.value
}
attr(default,"name")<- "default"
are<-function(object, class2,alls=T,...){
	if(is(object,'list')){zx<-sapply(object,FUN=Is,class2 = class2,...)}
	else{zx<-Is(object,class2= class2 ,...)}
	if(alls){zx<-all(zx==T)}
	zx
}
assert.is<-function (object, class2, text = "") {
    stopifnot(is.character(class2))
    if (missing(object)) 
        stop(paste(class2, "object missing in assert.is", text))
    stopifnot(length(class2) >= 1)
    if (!Is(object = object, class2 = class2)) {
        warning(paste("got ", class(object), sep = ""))
        message("got ", class(object), " when one of these classes was required:")
        print(class2)
        stop(text)
    }
}
assert.are <- function(object, class2, ...){
	assert.is(object, "list")
	if(!is.nothing(object)){
		for(obj in object)
			assert.is(object = obj, class2 = class2, ...)
	}
}
##--is\are---------
#  
is.prob <- function(P, tolerance = 1e-3, ...){
	boo <- is(P, "numeric") && all(P >= 0 - tolerance & P <= 1 + tolerance, ...)
	!is.na(boo) && boo
}

is.err <- function(object){all(is(object, "try-error"))}

#
is.ok<-function(object){all(object==T)}
is.valid<-function(object){
	k.valid<-function(object){!(length(object)==0) & !(is.na(object)) & !is(object,'try-error')}
	if(Is(object,c(d1class))){zo<-k.valid(object)}
	else if(is(object,'list')){zo<-sapply(object,k.valid)}
	zo}
is.vide<-function(object){length(object)==0}
is.error<-function(object){Is(object,c('try-error','NULL'))}
is.unique<-function(object){
	if(length(object)==0) {return(T)}
	assert.is(object,d1class)
	length(unique(object))==length(object)
}
#
#
any.mustbe.null<-function(x,y,...){
	zaux<-c(list(x,y),list(...))
	any(are.null(zaux,alls=FALSE)) && !are.null(zaux)
}
any.mustbe.nonnull<-function(x,y,...){
	zaux<-c(list(x,y),list(...))
	any(!are.null(zaux,alls=FALSE)) && !are.null(zaux)
}

##
are.prob <- function(object,alls=T, ...){
	if(is(object,'list')){zo<-sapply(object,is.prob,...)} else {zo<-is.prob(object,...)}
	if(alls){zo<-all(zo==T)}
	zo}
are.nothing <- function(object,alls=T){
	if(is(object,'list')){zo<- sapply(object, is.nothing)} else {zo<-is.nothing(object)}
	if(alls){zo<-all(zo==T)}
	zo}
are.null <- function(object,alls=T){
	if(is(object,'list')){zo<-sapply(object,is.vide)} else {zo<-is.vide(object)}
	if(alls){zo<-all(zo==T)}
	zo}
are.unique<-function(object,alls=T){
	if(is(object,'list')){zo<-sapply(object,is.unique)} else {zo<-is.unique(object)}
	if(alls){zo<-all(zo==T)}
	zo}
are.na<-function(object,alls=T){
	if(is(object,'list')){zo<-sapply(object,is.na)} else {zo<-is.na(object)}
	if(alls){zo<-all(zo==T)}
	zo}

are.error<-function(object,alls=T){
	if(is(object,'list')){zo<-sapply(object,is.error)} else{zo<-is.error(object)}
	if(alls){zo<-all(zo==T)}
	zo}
are.err<-function(object,alls=T){if(is(object,'list')){zo<-sapply(object,is.err)} else{zo<-is.err(object)}
	if(alls){zo<-all(zo==T)}
	zo}
are.ok<-function(object,alls=T){if(is(object,'list')){zo<-sapply(object,is.ok)} else{zo<-is.ok(object)}
	if(alls){zo<-all(zo==T)}
	zo}
are.valid<-function(object,alls=T){
	if(is(object,'list')){zo<-sapply(object,is.valid)} else{zo<-is.valid(object)}
	if(alls){zo<-all(zo==T)}
	zo}
#
is.named<-function(object){
 if( Is(object,d1class)){length(names(object))==length(object)}
else if(Is(object,d2class)){length(rownames(object))==nrow(object) && length(colnames(object))==ncol(object) }
  }	
are.named<-function(object,alls=T){
    if(is(object,'list')){zo<-sapply(object,is.named)}
	else {zo<-is.named(object)}
	if(alls){zo<-all(zo==T)}
	zo}

is.videna<-function(object){if(is.vide(object)){return(T)} else (is.na(object))}
are.nullna<-function(object,alls=T){
	if(is(object,'list')){zo<-sapply(object,is.videna)} else {zo<-is.videna(object)}
	if(alls){zo<-all(zo==T)}
	zo}


is.rowmatrix<-function(x){okm<-Is(x,d2class);if(!okm){return(FALSE)};oknr<-nrow(x)==1;if(okm&oknr){T} else {FALSE}}
is.colmatrix<-function(x){okm<-Is(x,d2class);if(!okm){return(FALSE)};oknc<-ncol(x)==1;if(okm&oknc){T} else {FALSE}}
##--------coversions---------



as.rowmatrix<-function(x){
	if (is.d1class(x)){zo<-matrix(x,1,length(x));colnames(zo)<-names(x);rownames(zo)<-'X1'}
	else if (Is(x,d2class) && nrow(x)==1){zo<-x}
	else{zo<-x}
	zo
}
as.colmatrix<-function(x){
	if (is.d1class(x)){zo<-matrix(x,length(x),1);rownames(zo)<-names(x);colnames(zo)<-'I1'}
	else if (Is(x,d2class) && ncol(x)==1){zo<-x}
	else{zo<-x}
	zo
}

as.num<-function(x){
	if(is.rowmatrix(x)){nm<-colnames(x)}
	if(is.colmatrix(x)){nm<-rownames(x)}
	zo<-as(x,'numeric');names(zo)<-nm
	zo
}

factor2numeric<-function(x){
    if(is(x,'numeric')){return(x)}
    else if(is(x,'factor')){as.numeric(levels(x)[x])}}
factor2character<-function(x){
    if(is(x,'character')){return(x)}
    else if(is(x,'factor')){as.character(levels(x)[x])}}

est2list<-function(x){
  csok<-c(d1class,d2class,'scalar','NULL','function','list')
  if (class(x)%in%csok){return(x)}
  if (class(x)=='list'){zx<-lapply(x,est2list);return(zx)}
  #isS4(object)
  slnm<-slotNames(x)
 
 if(length(slnm)==0){return(x)}
 
   z<-lapply(1:length(slnm),function(i){
	z1<-eval(parse(text=paste('x@',slnm[i],sep='')))	
	if(class(z1)%in%csok){return(z1)}
	else {z1<-est2list(x=z1)}
	z1})
   names(z)<-slnm
   assert.is(z,'list')
   class(z)<-c('list',paste(class(x),'list',sep='.'))
   z
}
list2est<-function(x,n.object){
    csok<-c(d1class,d2class,'scalar','NULL','function')#'Numeric','Scalar','Matrix',
    cz<-NULL
	if(is(x,'list') & length(class(x))>1) {
		ind<-grep(x=class(x),pattern='.list',fixed=T)
		if(length(ind)==0){break()}
		obj<-class(x)[ind]
		n.object<-string2vect(obj,sep='.list')[1]}
		
	if(!is(x,'list') || is.vide(n.object)){return(x)}
	nms<-slotNames(n.object)
	assert.is(x,'list')
	if(!all(nms%in%names(x))){if(echo==TRUE){message('could not convert list to ',n.object)};return(x)}
	nx<-new(n.object)
	for (i in 1:length(nms)){
		eval(parse(text=paste('cz<-class(nx@',nms[i],')',sep='')))#XXX|:no visible binding for global variable ÔczÕ
		if(cz%in%csok){eval(parse(text=paste('nzi<-list2est(x$',nms[i],',n.object=cz)',sep='')))}
		else{eval(parse(text=paste('nzi<-x$',nms[i],sep='')))}
		eval(parse(text=paste('nx@',nms[i],'<-nzi',sep='')))
  }
  validObject(nx)
  nx
}
#
data.frame2list<-function(x){#as.list(x)
    assert.is(x,c('data.frame','matrix',d2class))
    rnm<-rownames(x)
    zo<-lapply(1:ncol(x),FUN=function(j){z<-x[,j];names(z)<-rnm;z})
    names(zo)<-colnames(x)
    zo
}
list2data.frame<-function(x){as.data.frame(x,stringsAsFactors=FALSE)}

matrix2list<-function(x,row=T,...){
	if(Is(x,d1class)){x<-as.rowmatrix(x)}
	x<-MakeNames(x,...);cnm<-colnames(x);rnm<-rownames(x)
	if(!row||is(x,'data.frame')){zo<-data.frame2list(x)}
	else if(row &&!is(x,'data.frame')){
		zo<-lapply(1:nrow(x),FUN=function(j){z<-x[j,];names(z)<-cnm;z})
		names(zo)<-rownames(x)
	}
	zo
}
list2matrix <- function(x,row=T,...){
	if(Is(x,c(d1class,d2class))){return(x)}
	if(isS4(x)){x<-est2list(x)}
	assert.is(x,'list')
	zo<-c();nms<-c()
	
	for(i in 1:length(x)){z<-x[[i]]
		if(!Is(z,c(d1class,d2class))){next()}
		if(row){zo<-rbindx(zo,x[[i]]);nms<-c(nms,names(x)[i])}
		else{zo<-cbindx(zo,x[[i]]);nms<-c(nms,names(x)[i])}
	}
	if(row){rownames(zo)<-nms} else {colnames(zo)<-nms}

	  zo
}
#comodidad---------
eval.txt<-function(txt){z<-NULL;eval(parse(text=paste('z<-',txt,sep='')));z}#XXX|:no visible binding for global variable ÔzÕ
npaste<-pastex<-function(...){paste(...,sep='')}
nmin<-minx<-function(x,...){min(x,...,na.rm=TRUE)}
nmax<-maxx<-function(x,...){max(x,...,na.rm=TRUE)}
nmean<-meanx<-function(x,...){mean(x,...,na.rm=TRUE)}
nmedian<-medianx<-function(x,...){median(x,...,na.rm=TRUE)}
ndata.frame<-data.framex<-function(...){data.frame(...,stringsAsFactors=FALSE)}
nsd<-sdx<-function(x,...){sd(x,...,na.rm=TRUE)}
nmapply<-mapplyx<-function(...){mapply(...,SIMPLIFY=FALSE)}
#matrix-----
undefAsNA<-function(x){
if(is(x,'numeric')){ind<-which(!is.finite(x));if(length(ind)>0){x[ind]<-NA} else {x}}
else if(is(x,'matrix')){x<-apply(x,2,undefAsNA)}
x
}
k.removeRC.from.matrix<-function(x,opt='NA',indx=FALSE){
	
	opt<-opt[opt%in%c('NA','na','eq','eqs','EQ','EQS','nfinit','Inf','inf','NaN','nan')]
	if(length(opt)==0){stop('bad opt')}
	indr<-indc<-c()
	for(i in 1:length(opt)){
		z<-k.removeRC.from.matrix(x=x,opt=opt[i],indx=T)
		x<-z$matrix
		indr<-c(indr,z$indNA.rows)
		indc<-c(indc,z$indNA.cols)
	}
	if(indx){list(matrix=x,indNA.rows=indr,indNA.cols=indc)}
	else{x}
}

removeRC.from.matrix<-function(x,opt='NA',indx=FALSE){ox<-x;if(is.vide(x)){message('empty matrix');return(NULL)}
	if (Is(x,d1class)){x<-as.rowmatrix(x)}
	    all.eq<-function(x){if(length(table(x))==1){1} else{0}
	}
	all.na<-function(x){
		if(length(which(is.na(x)))==length(x)){1} else{0}
	}
	all.nonfinit<-function(x){
		if(length(which(!is.finite(x)))==length(x)){1} else{0}}
	zr.fun<-function(x,fun){
		assert.is(x,d2class)
		zr<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){fun(x[i,])})
		names(zr)<-rownames(zr);zr
	}
	zc.fun<-function(x,fun){assert.is(x,d2class)
		zr<-vapply(1:ncol(x),FUN.VALUE=numeric(1),FUN=function(i){fun(x[,i])})
		names(zr)<-colnames(zr);zr
	}
	rem.row<-function(x,zr){
		if(1 %in% zr){x<-x[-which(zr==1),,drop=FALSE]}
		x}
	rem.col<-function(x,zc){
		if(1 %in% zc){	x<-x[,-which(zc==1),drop=FALSE]}
		x
	}

	if(opt%in%c('NA','na'))	{zr<-zr.fun(x=x,fun=all.na);zc<-zc.fun(x=x,fun=all.na)}
	else if(opt%in%c('eq','eqs','EQ','EQS')){zr<-zr.fun(x=x,fun=all.eq);zc<-zc.fun(x=x,fun=all.eq)}
	else if(opt%in%c('nfinit','Inf','inf','NaN','nan')){zr<-zr.fun(x=x,fun=all.nonfinit);zc<-zc.fun(x=x,fun=all.nonfinit)}
	else{stop('bad opt')}
	
x<-rem.row(x=x,zr=zr);indr<-which(zr==1)
x<-rem.col(x=x,zc=zc);indc<-which(zc==1)
if(Is(ox,d1class)){nms<-colnames(x);x<-as.numeric(x);names(x)<-nms}
	if(echo){
message(paste('removing',length(which(zr==1)),'rows and',length(which(zc==1)),'columns\n') )}
if(indx){list(matrix=x,indNA.rows=indr,indNA.cols=indc)}
else{x}
}

removeNA.from.matrix<-function(x,indx=FALSE){removeRC.from.matrix(x=x,opt='nfinit',indx=indx)}#x<-undefAsNA(x);
removeEQ.from.matrix<-function(x,indx=FALSE){removeRC.from.matrix(x=x,opt='eq',indx=indx)}
removeNF.from.matrix<-function(x,indx=FALSE){removeRC.from.matrix(x=x,opt='nfinit',indx=indx)}

#
nsort<-function(x,ind=FALSE,...){
	ok<-is.rowmatrix(x)|is.colmatrix(x)
	if(ok){x<-as.num(x)}
	zo<-sort(x=x, index.return=TRUE,...)
	if(ind){return(zo$ix)}
	else{return(zo$x)}
}

k.indSort<-function(x,ref,inter=FALSE,ind=T){y<-ref;ox<-x;oy<-y
	if(length(ref)==0){stop('indSort: ref is NULL')}
	if(length(x)==0){stop('indSort: x is NULL')}

	x[which(is.na(x))]<-'NA'
	y[which(is.na(y))]<-'NA'
	x <- make.unique(as.character(x))
	y <- make.unique(as.character(y))
	
	xx<-c(1:length(x));yy<-c(1:length(y))
	names(xx)<-names(x)<-x;names(yy)<-names(y)<-y  
	
	finty<-intersect(y,x)
	if(length(finty)==0){stop('input vectors do not have common names')}
		#if(echo){message('k.indSort: x and ref have nothing in common, returning NULL')};return(NULL)}
	if(inter) {
		fintx<-finty
		nms<-finty}
	else{
		indx<-finty
		indnx<-which(!x%in%y)
		nn<-length(x)+length(y)-length(finty)
		w<-rep(NA,nn)
		
		names(w)[1:length(y)]<-y
		#if(length(indnx)>0){}
		indw<-c((length(y)+1):(length(y)+length(indnx)))
		indw<-indw[which(indw>=0 & indw<=(length(y)+length(indnx)) & indw>=indw[1])]
		stopifnot(length(indw)==length(indnx))
		names(w)[indw]<-x[indnx]
		nms<-names(w)
		
		
		ny<-w;ny[1:length(y)]<-yy
		nx<-w;nx[indx]<-xx[indx];nx[indw]<-xx[indnx]
		
		kny<-w;kny[1:length(y)]<-y
		knx<-w;knx[indx]<-indx;knx[indw]<-x[indnx]
		
		names(kny)<-names(knx)<-NULL
		names(ny)<-names(nx)<-NULL
		fintx<-knx
		finty<-kny
	}
	  
	
	  indsortedx<-xx[fintx];names(indsortedx)<-nms
	  indsortedy<-yy[finty];names(indsortedy)<-nms
	 
	vfintx<-fintx[which(!is.na(fintx))]
	vfinty<-finty[which(!is.na(finty))]
	if(length(vfintx)>0 && length(vfinty)>0){
		ok0<-all(x[vfintx]==vfintx)&&all(y[vfinty]==vfinty)
		ok3<-all(y[indsortedy][which(!is.na(y[indsortedy]))]==y[indsortedy][which(!is.na(y[indsortedy]))])
	}
	if(length(vfintx)>0){
		ok1<-all(x[indsortedx][which(!is.na(x[indsortedx]))]==vfintx)
	}
	if(length(vfinty)>0){
		ok2<-all(y[indsortedy][which(!is.na(y[indsortedy]))]==vfinty);ok1<-ok0<-T	
	}
	if(length(vfintx)>0 && length(vfinty)==0){ok2<-ok0<-ok3<-T}
	if(length(vfintx)==0 && length(vfinty)>0){ok1<-ok0<-ok3<-T}

	stopifnot(ok0&&ok1&&ok2&&ok3)
	  if(ind){
		list(x=indsortedx,ref=indsortedy,names=nms)}
	else{
		nx<-ox[indsortedx];names(nx)<-nms
		ny<-oy[indsortedy];names(ny)<-nms
		list(x=nx,ref=ny,names=nms)}
}

indSort<-function(x,ref=NULL,ind=T,rm.na=FALSE,inter=FALSE){
	
	y<-ref;xx<-x;yy<-y
	if(length(ref)==0){stop('indSort: ref is NULL')}
	if(length(x)==0){stop('indSort: x is NULL')}

	fok<-intersect(y,x)
	if(length(fok)==0&&inter){warning('indSort: x and ref have nothing in common, returning NULL');return(NULL)}
	z<-k.indSort(x=x,ref=ref,ind=T)
	indsorted<-z$x;indy<-z$ref
        
	nindx<-NA*c(1:(length(y)+length(which(!x%in%fok))))
	nindx[which(y%in%fok)]<-indsorted 
	
	indsorted<-nindx

	stopifnot(identical(x[indsorted[is.finite(indsorted)]],fok))
	
	if(rm.na | inter){indsorted<-indsorted[is.finite(indsorted)]}
	
	if(ind & !inter){zo<-list(x=indsorted,ref=1:length(y))}
	else if(!ind & !inter){zo<-list(x=xx[indsorted],ref=yy)}
	else if(ind & inter){zo<-list(x=indsorted,ref=indy)}
	else if(!ind & inter){zo<-list(x=xx[indsorted],ref=yy[indy])}
	
zo	
}
indSortAsY<-function(x,y,inter=T,ind=T){

	z<-k.indSort(x=x,ref=y,ind=ind,inter=inter)
	if(length(z)==0){return(NULL)}
	names(z)[which(names(z)%in%'ref')]<-'y'
	names(x)<-names(y)<-NULL
	
	if(inter&&ind){ok<-identical(x[z$x],y[z$y])}
	else if(inter&&!ind){ok<-identical(z$x,z$y)}
	else if(!inter){ok<-identical(names(z$x),names(z$y))}
	stopifnot(ok)

	z
}
#
make.labels<-function(n,nmvar=c('X'),n.ini=1){n.ini<-n.ini[1]
    if(is.vide(n)){return(NULL)}
    stopifnot(n>0)
    paste(nmvar[1],c(1:n)+n.ini-1,sep='')
}



k.bind<-function(x,y,ind=T,inter=FALSE){
  if(length(x)==0){return(y)}
  if(length(y)==0){return(x)}
      assert.are(list(x,y),d1class)
	ox<-x;oy<-y
	x<-make.unique(as.character(x));y<-make.unique(as.character(y))
	fx<-intersect(x,y)
	if(inter&ind){
		z<-k.indSort(x=y,ref=x,ind=T,inter=T)
		stopifnot(identical(x[z$ref],y[z$x]))
		xx<-z$ref
		yy<-z$x
		xy<-x[z$ref]
	}
	else if(!inter&!ind){
		xx<-c(fx,setdiff(x,fx),rep(NA,length(setdiff(y,fx))))
		yy<-c(fx,rep(NA,length(setdiff(x,fx))),setdiff(y,fx))
		xy<-c(fx,setdiff(x,fx),setdiff(y,fx))	
		}
	else if(!inter&ind){
		z<-k.indSort(x=y,ref=x,ind=T,inter=T)
		xx<-c(z$ref,which(x%in%setdiff(x,fx)),rep(NA,length(setdiff(y,fx))))
		yy<-c(z$x,rep(NA,length(setdiff(x,fx))),which(y%in%setdiff(y,fx)))
		xy<-c(x[z$ref],x[which(x%in%setdiff(x,fx))],y[which(y%in%setdiff(y,fx))])
	}

	stopifnot(length(xx)==length(yy)&& length(xx)==length(xy))
	
list(x=xx,y=yy,all=xy)
	
}
k.MakeNames<-function(x,col=T,nmvar='I',force=FALSE,unique=T,n0=1){
		if(is.vide(x)){return(x)}
	k.nms<-function(nm,nn){
		if(is.vide(nm) | all(is.videna(nm)) | force==T ){
			nm<-make.labels(n=nn,nmvar=nmvar[1],n.ini=n0[1])}
		if(any(is.na(nm))){ind<-which(is.na(nm));nm[ind]<-paste(nmvar[1],ind,sep='')}
		if(unique==T){nm<-make.unique(nm)}
		nm
			
	}
	if(Is(x,c(d1class,'list'))){
		nm<-k.nms(nm=names(x),nn=length(x))
		#if(is.vide(nm) | all(is.videna(nm)) | force==T){
		#	nm<-make.labels(n=length(x),nmvar=nmvar[1],n.ini=n0[1])}
		#if(unique==T){nm<-make.unique(nm)}
		#if(any(is.na(nm))){ind<-which(is.na(nm));nm[ind]<-paste(nmvar[1],ind,sep='')}
		names(x)<-nm;return(x)}
	
	assert.is(x,d2class)
	nn<-nsize(x)
	
	if(col){nm<-colnames(x);nn<-nn[2]}
	else{nm<-rownames(x);nn<-nn[1]}

	nm<-k.nms(nm=nm,nn=nn)
	#if(all(is.vide(nm) |is.videna(nm)) | force==T) {nm<-make.labels(n=nn,nmvar=nmvar[1],n.ini=n0[1])}
	#if(unique==T){nm<-make.unique(nm)}
	
	if(col){colnames(x)<-nm}
	else{rownames(x)<-nm}
	x
}
MakeColNames<-function(x,...){k.MakeNames(x,col=T,...)}
MakeRowNames<-function(x,...){k.MakeNames(x,col=FALSE,...)}
MakeNames<-function(x,nmvar=c('X','I'),...){
	if(Is(x,c(d1class,'list'))){x<-k.MakeNames(x=x,nmvar=nmvar[1],...)}
	else if(Is(x,d2class)){
		x<-k.MakeNames(x,col=T,nmvar=nmvar[pmax(2,length(nmvar))],...)
		x<-k.MakeNames(x,col=FALSE,nmvar=nmvar[1],...)
	}
	x
}
k.add.NA.col.dataframe<-function(x,ind){
	if(Is(x,d1class)){zo<-x[ind];return(zo)}
	assert.is(x,d2class)
	if(all(is.finite(ind)==T)){zo<-x[,ind,drop=FALSE];return(zo)}
	if(is(x,'matrix')){zo<-x[,ind,drop=FALSE];return(zo)}
	assert.is(x,'data.frame')
	names(ind)<-1:length(ind)
	indna<-which(!is.finite(ind))
	z1<-cbind(x[,ind[is.finite(ind)],drop=FALSE],matrix(NA,nrow(x),length(indna)))
	colnames(z1)<-c(colnames(x)[ind[is.finite(ind)]],paste('Unk',1:length(indna),sep=''))
	nms<-c(names(ind[is.finite(ind)]),names(ind[indna]))
	z<-k.indSort(x=nms,ref=names(ind),ind=T,inter=T)
	z1[,z$x,drop=FALSE]
}

k.nrbind<-k.rbindx<-function(x,y,inter=FALSE){
  if(length(x)==0){return(y)}
  if(length(y)==0){return(x)}
  if(identical(pnames(x),pnames(y))){zo<-rbind(x,y);return(zo)}
	if(Is(x,d1class)){x<-as.rowmatrix(x)}
	if(Is(y,d1class)){y<-as.rowmatrix(y)}

	x<-MakeColNames(x,force=FALSE,unique=T)
	y<-MakeColNames(y,force=FALSE,unique=T)
	
	if(identical(colnames(x),colnames(y))){return(rbind(x,y))}
	z<-k.bind(x=colnames(x),y=colnames(y),ind=T,inter=inter)

	z1<-k.add.NA.col.dataframe(x=x,ind=z$x)
	rownames(z1)<-rownames(x)
	z2<-k.add.NA.col.dataframe(x=y,ind=z$y)
	rownames(z2)<-rownames(y)
	colnames(z1)<-colnames(z2)<-z$all

	rbind(z1,z2)
  
}
k.ncbind<-k.cbindx<-function(x,y,inter=FALSE){
  if(length(x)==0){return(y)}
  if(length(y)==0){return(x)}
  if(identical(fnames(x),fnames(y))){zo<-cbind(x,y);return(zo)}
  
	if(Is(x,d1class)){x<-as.colmatrix(x)}
	if(Is(y,d1class)){y<-as.colmatrix(y)}
	
	x<-MakeRowNames(x,force=FALSE,unique=T)
	y<-MakeRowNames(y,force=FALSE,unique=T)

	if(identical(rownames(x),rownames(y))){return(cbind(x,y))}
	z<-k.bind(x=rownames(x),y=rownames(y),ind=T,inter=inter)

	z1<-x[z$x,,drop=FALSE];colnames(z1)<-colnames(x)
	z2<-y[z$y,,drop=FALSE];colnames(z2)<-colnames(y)
	zo<-cbind(z1,z2)
	rownames(zo)<-z$all
	zo
}
ncbind<-cbindx<-function(x,y,inter=FALSE){
  if(length(x)==0){return(y)}
  if(length(y)==0){return(x)}
	if(are(list(x,y),c(d1class,d2class))){k.ncbind(x=x,y=y,inter=inter)}
	else if(are(list(x,y),'list')){
	   x<-MakeNames(x=x,force=FALSE,unique=T);y<-MakeNames(x=y,force=FALSE,unique=T)
		z<-k.indSort(x=names(y),ref=names(x),ind=T,inter=T)
		nx<-x[z$ref];ny<-y[z$x]
		stopifnot(identical(names(nx),names(ny)))
		zo<-mapply(x=nx,FUN=cbindx,y=ny,inter=inter,SIMPLIFY=FALSE)
		return(zo)}
}
nrbind<-rbindx<-function(x,y,inter=FALSE){
    if(length(x)==0){return(y)}
  if(length(y)==0){return(x)}
  if(are(list(x,y),c(d1class,d2class))){k.nrbind(x=x,y=y,inter=inter)}
  else if(are(list(x,y),'list')){
    x<-MakeNames(x=x,force=FALSE,unique=T);y<-MakeNames(x=y,force=FALSE,unique=T)
    z<-k.indSort(x=names(y),ref=names(x),ind=T,inter=T)
		nx<-x[z$ref];ny<-y[z$x]
		stopifnot(identical(names(nx),names(ny)))
		zo<-mapply(x=nx,FUN=rbindx,y=ny,inter=inter,SIMPLIFY=FALSE)
		return(zo)}
}
###======================================================
same.length<-same.size<-same.dim<-function(x,y,exact=FALSE){
	if(exact){identical(nsize(x),nsize(y))}
	else{setequal(nsize(x),nsize(y))}
}
SameNames<-function(x,y,exact=FALSE){
	if(are(list(x,y),d1class)||are(list(x,y),'list')){mx<-names(x);my<-names(y)}
	if(are(list(x,y),d2class)){mx<-c(rownames(x),colnames(x));my<-c(rownames(y),colnames(y))}
	if(exact){identical(mx,my)}
	else{setequal(mx,my)}
}
same.content<-function(x,y=NULL,tol=1e-10){  
    k.unit.same.vectors<-function(x,y,tol){names(x)<-names(y)<-NULL
	if(!same.length(x=x,y=y,exact=T)){return(FALSE)}
	if(!SameNames(x=x,y=y,exact=T)){warning('not same names')}
	#if(!same.name(x=x,y=y,exact=T)){return(FALSE)}
	if(are(list(x,y),c("numeric"))){z<-abs(x-y);all(z[!is.na(z)]<=tol)}
	else if(are(list(x,y),c("character","logical"))){identical(x,y)}
	else{stop('not yet implemented')}}
if((is.vide(x)&!is.vide(y)) || (is.vide(y)&!is.vide(x))){return(FALSE)}
else if(are.null(list(x,y))){return(T)}

if(are(list(x,y),d1class) ){
	zo<-k.unit.same.vectors(x=x,y=y,tol=tol)}
else if (are(list(x,y),d2class) ){
	zo<-apply(x,2,FUN=k.unit.same.vectors,y=y,tol=tol)}
else if (are(list(x,y),'list') ){
	zo<-unlist(lapply(x,FUN=k.unit.same.vectors,y=y,tol=tol))}
	all(zo==T)
}


###--functions------------------------------------------------------------------
Is.function<-function(x){i=0;n=5
	if(!Is(x,c('function','character','standardGeneric'))){return(FALSE)}
	if(is.function(x)){return(TRUE)}
	while(is(x,'character')&&i<n){
		z<-try(eval(parse(text=x)),silent=T);x<-z
		if(is.error(x) | is.function(x) | is(x,'standardGeneric')){break()}
		i<-i+1}
	is.function(x)}
As.function<-function(x){i=0;n=5
	assert.is(x,c('function','character','standardGeneric'))
	if(is(x,'function')){return(x)}
	if(!Is.function(x)){stop('As.function: bad x')}
	while(is(x,'character')&&i<n){
		z<-eval(parse(text=x));x<-z
		if(is.error(x) || is.function(x)|| is(x,'standardGeneric')){break()}
		i<-i+1}
	return(x)}

setAs(from='character',to='function',function(from){
    z<-try(eval(parse(text=from)),silent=FALSE)
    if(is.error(z)){message('object ',from,' is not a function');return(NULL)}
    z})

getName<-function(x){
	z2<-as.list(match.call(expand.dots=TRUE)[-1])[[1]]
	z2<-as.character(z2)
        if(Is.function(x)){x<-As.function(x=x)}
	zo<-attr(x,'name')
	if(!is.vide(zo)){return(zo)}
        if(is(x,'character')){return(x)}
	else{return(z2)}
	# else if(is(x,'list')){zo<-lapply(x,getName)}
        
}
setGeneric("setName", function(x,...) standardGeneric("setName"))
#setName<-function(x,value){}
setMethod("setName", signature(x='character' ),function(x,value=NULL){
	if(is.vide(value)){value<-x}
	#if(!Is.function(x)){stop('could not name x')}
	z<-try(eval(parse(text=paste('attr(',x,',"name")<-"',value,'"',sep='')),envir =.GlobalEnv),silent=T)
	if(is.err(z)){stop('could not name x')}})








NameObjectsInRpackage<-function(Rpack='stats',mode='function',pattern='.test'){
  za<-ls.str(paste("package:",Rpack,sep=''), mode = mode,pattern=pattern)
  za<-unlist(as.list(za))
  for(i in 1:length(za)){
	z<-try(eval(parse(text=paste('attr(',za[i],',"name")<-"',za[i],'"',sep='')),envir =.GlobalEnv),silent=T)
	if(is.err(z)){next()}
	}}
   
any2char<-function(x){
	if(!is.function(x)){as.character(x)}
	cnm<-as.character(as.list(match.call(definition =x))[[2]])}
getName<-function(x){zo<-NULL;ox<-x
	nm<-attr(x, "name");if(length(nm)>0){return(nm)}
	if(is.character(x)){nx<-try(As.function(x),silent=T);if(is.err(nx)){return(x)} else{x<-nx}}
	if(Is.function(x)){
		cnm<-as.character(as.list(match.call(definition =x))[[2]])
		
		ob<-ls.str(, mode = "function",envir=environment(x))
		ob<-unlist(as.list(ob))
		if(is.character(ox) || cnm%in%ob){if(ox%in%ob){return(ox)} else if(ox%in%ob){return(cnm)}}
		
		fob<-mget(ob, envir = environment(x), mode = "function",ifnotfound=NA)
		zok<-vapply(1:length(ob),FUN.VALUE=logical(1),FUN=function(i){identical(functionBody(x),functionBody(fob[[i]]))})
		names(zok)<-names(fob)
		if(length(zok)==0 || length(which(zok==T))==0){if(is.character(ox)){return(ox)} else {return(cnm)}}
		#else if(length(zok)>1){if(is.character(ox)){return(ox)} else {return(cnm)}}
		zo<-names(zok)[zok==T];return(zo)}
	zo
}
get.args<-function(fun,max.size=c(1,3),n=1,wn=-1,...){
  z1<-try(match.call(definition =fun,call=sys.call(which = wn)),silent=T)
  z2<-try(formals(fun = fun),silent=T)
  z2<-z2[!names(z2)%in%"..."]
  zz<-try(uniquex(c(as.list(z1[2:length(z1)]),as.list(z2),list(...))),silent=T)
 
  infox<-lapply(1:length(zz),FUN=function(i){
	#if(is(zz[[i]],"name")){zz[[i]]<-as.character(zz[[i]])}
	zx<-try(eval(zz[[i]],envir=parent.frame(n=n)),silent=T)
	if(is.function(zx)){zx<-getName(zx)}
	else if(is.err(zx)){zx<-names(zz[i])}
	if(!Is(zx,c(d1class,d2class))){zx<-class(zx)}
	zx});names(infox)<-names(zz)
sz<-sapply(infox,length)
  c(infox[sz>=max.size[1] & sz<max.size[2]],fun.name=getName(fun))
}




##============xprnSet=======================

setClass("XprnSetPair", representation("xprnSetPair"))
setValidity("XprnSetPair", function(object){
	xx<-exprs(object@x);yy<-exprs(object@y)
	xx<-as.numeric(xx);yy<-as.numeric(yy)
	all(xx[is.finite(xx)] >= 0, na.rm = TRUE)
	all(yy[is.finite(yy)] >= 0,  na.rm = TRUE)
})

setClassUnion(name = "XprnSetObject", members = c("XprnSet", "XprnSetPair"))
setClassUnion(name = "prnSet", members = c("XprnSet", "xprnSet"))
setClassUnion(name = "prnSetPair", members = c("XprnSetPair", "xprnSetPair"))
#-------

setClassUnion("missingOrNULL", c("missing", "NULL"))
setClassUnion("matrixORnumeric", c("matrix","numeric"))
setClassUnion("MatrixORNumeric", c("Matrix","Numeric"))
setClassUnion("all.matrixORnumeric", c("Matrix","matrix","numeric","Numeric"))#"character","matrix","missing", "NULL",
setClassUnion("pfdata", c("data.frame","AnnotatedDataFrame"))#"character","matrix","missing", "NULL",
setClassUnion("missingOrNULLORpfdata", c("data.frame","AnnotatedDataFrame","missing", "NULL"))
setClassUnion("all.pfdata", c("AnnotatedDataFrame",d1class,d2class,"missing", "NULL"))



setAs(from = "matrixORnumeric", to = "Matrix", function(from){Matrix(from)})
setAs(from = "numeric", to = "Matrix", function(from){Numeric(from)})
as.Matrix<-function(from){
	if(class(from)%in%c('Matrix','Numeric')){from}
	if(Is(from,c("matrix","numeric"))){Matrix(from)}	
}
as.matrices<-function(from){
	if(class(from)%in%c('Matrix')){as(from,'matrix')}
	else if(class(from)%in%c('Numeric')){as(from,'numeric')}
	else{ as.matrix(from)}	
}
#-------
prnSet<-c("XprnSet", "xprnSet")
prnSetPair<-c("XprnSetPair", "xprnSetPair")
is.xprnSet<-function(x){class(x)%in%c("xprnSet")}
is.XprnSet<-function(x){class(x)%in%c("XprnSet")}
is.prnSet<-function(x){is.xprnSet(x)||is.XprnSet(x)}
is.xprnSetPair<-function(x){class(x)%in%c("xprnSetPair")}
is.XprnSetPair<-function(x){class(x)%in%c("XprnSetPair")}
is.prnSetPair<-function(x){is.xprnSetPair(x)||is.XprnSetPair(x)}
is.Matrix<-function(x){Is(x,c('matrix','numeric'))&&class(x)==c('Matrix','Numeric')}
nis.Matrix<-function(x){Is(x,c('matrix','numeric'))&&!class(x)%in%c('Matrix','Numeric')}
is.ExpressionSet<-function(x){is(x,'ExpressionSet')&&!Is(x,c(prnSet,prnSetPair))}
is.matrixORnumeric<-function(x){Is(x,c("matrix","numeric"))}
is.MatrixORNumeric<-function(x){is.Matrix(x)}
#--------------------------------------
setMethod("annotation", signature(object = "prnSetPair"), function(object){
	paste(annotation(object@x), annotation(object@y), sep = " vs. ")
})


setReplaceMethod("annotation", signature(object="prnSet"),function(object, value) {
	z<-as(object, "ExpressionSet")
	 annotation(z)<-value
	 if(is(z,"XprnSet")){zo<-as.XprnSet(z)}
	 else if(is(z,"xprnSet")){zo<-as.xprnSet(z)}
zo})
setReplaceMethod("annotation", signature(object="prnSetPair"),function(object, value) {
	if(nchar(value)>0){
	annotation(object@x)<-paste(annotation(object@x),' -',value,'-',sep='')
	annotation(object@y)<-paste(annotation(object@y),' -',value,'-',sep='')}
	object})

#
setMethod("exprs", signature(object = "prnSetPair"), function(object){list(x=exprs(object@x),y=exprs(object@y))})
#
setMethod("is.paired", signature(x = "prnSetPair",y="missingOrNULL"), function(x,y){
        ax<-annotation(x)
        z0<-grep(pattern=paste('-pairedTRUE',sep=''),x=ax,fixed =TRUE)
        z1<-identical(colnames(exprs(x@x)),colnames(exprs(x@y)))
	
        if(length(z0)>0||z1==TRUE){return(TRUE)}
        else if(length(z0)==0||z1==FALSE){return(FALSE)}})
setMethod("is.paired", signature(x = "matrixORnumeric",y="missingOrNULL"), function(x,y){NULL})
setMethod("is.paired", signature(x = "prnSet",y="missingOrNULL"), function(x,y){NULL})
setMethod("is.paired", signature(x = "matrixORnumeric",y="matrixORnumeric"), function(x,y){
	
	if(is(x,'numeric')){x<-as.rowmatrix(x)}
	if(is(y,'numeric')){y<-as.rowmatrix(y)}
	assert.are(list(x,y),c('matrix'))
	x<-MakeNames(x)
	y<-MakeNames(y)
	if(identical(rownames(x),rownames(y)) && identical(colnames(x),colnames(y))){return(TRUE)}
		else{return(FALSE)}
		
	})
setMethod("is.paired", signature(x = "prnSet",y="prnSet"), function(x,y){
	ok1<-is.xprnSet(x)&&is.xprnSet(y) || is.XprnSet(x)&&is.XprnSet(y)
	ok2<-is.paired(x=exprs(x@x),y=exprs(x@y))
	if(!ok1){message('inconsistent');return(FALSE)}
		ok1&&ok2
	})
setMethod("is.positive", signature(x = "prnSet"), function(x){is.positive(exprs(x))})
setMethod("is.positive", signature(x = "prnSetPair"), function(x){is.positive(exprs(x@x))&&is.positive(exprs(x@y))})
#
setMethod("pData", signature(object = "prnSet"), function(object){
  pData(as(object, "ExpressionSet"))
})
setMethod("fData", signature(object = "prnSet"), function(object){
  fData(as(object, "ExpressionSet"))
})
setMethod("fData", signature(object = "prnSetPair"), function(object){
  zx<-fData(as(object@x, "ExpressionSet"))
  zy<-fData(as(object@y, "ExpressionSet"))
  stopifnot(identical(rownames(zx),rownames(zy)))
  zz<-cbindx(zx,zy)
  
  #z1<-as.data.frame(zz,stringsAsFactors=FALSE);rownames(z1)<-rownames(zx)
zz})
setMethod("pData", signature(object = "prnSetPair"), function(object){
  zx<-pData(as(object@x, "ExpressionSet"))
  zy<-pData(as(object@y, "ExpressionSet"))
  list(x=zx,y=zy)
})
setReplaceMethod("pData", signature(object="prnSet"),function(object, value) {
	z<-as(object, "ExpressionSet")
	 pData(z)<-value
	 if(is(object,"XprnSet")){zo<-as.XprnSetx(z)}
	 else if(is(object,"xprnSet")){zo<-as.xprnSetx(z)}
	 zo})
setReplaceMethod("fData", signature(object="prnSet"),function(object, value) {
	z<-as(object, "ExpressionSet")
	 fData(z)<-value
	 if(is(object,"XprnSet")){zo<-as.XprnSetx(z)}
	 else if(is(object,"xprnSet")){zo<-as.xprnSetx(z)}
	 zo})
#---
setMethod("colnames", signature(x='prnSet' ),function(x){colnames(exprs(x))})
setMethod('ncol', signature(x = "prnSet" ), function(x){ncol(exprs(x))})
setMethod("rownames", signature(x='prnSet' ),function(x){rownames(exprs(x))})
setMethod('nrow', signature(x = "prnSet" ), function(x){nrow(exprs(x))})
setMethod("colnames", signature(x='prnSetPair' ),function(x){list(x=colnames(exprs(x@x)),y=colnames(exprs(x@y)))})
setMethod('ncol', signature(x = "prnSetPair" ), function(x){list(x=ncol(exprs(x@x)),y=ncol(exprs(x@y)))})
setMethod("rownames", signature(x='prnSetPair' ),function(x){rownames(exprs(x@x))})
setMethod('nrow', signature(x = "prnSetPair" ), function(x){nrow(exprs(x@x))})
setMethod("fnames", signature(object = "prnSet"), function(object){featureNames(as(object,"ExpressionSet"))})
setMethod("pnames", signature(object = "prnSet"), function(object){colnames(exprs(object))})
setMethod("pnames", signature(object = "prnSetPair"), function(object){list(x=colnames(exprs(object@x)),y=colnames(exprs(object@y)))})
setMethod("fnames", signature(object = "prnSetPair"), function(object){featureNames(as(object@x,"ExpressionSet"))})

setReplaceMethod("fnames", signature(object="prnSet"),function(object, value) {
	z<-as(object, "ExpressionSet")
	 featureNames(z)<-value
	 if(is(object,"XprnSet")){zo<-as.XprnSetx(z)}
	 else if(is(object,"xprnSet")){zo<-as.xprnSetx(z)}
	 zo})
setReplaceMethod("pnames", signature(object="prnSet"),function(object, value) {
	z<-exprs(object)
	colnames(z)<-value
	exprs(object)<-z
	object})
setReplaceMethod("fnames", signature(object="prnSetPair"),function(object, value) {
	zx<-as(object@x, "ExpressionSet");zy<-as(object@y, "ExpressionSet")
	 featureNames(zx)<-featureNames(zy)<-value
	 if(are(list(object@x,object@y),"XprnSet")){zx<-as.XprnSetx(object@x);zy<-as.XprnSetx(object@y)}
	 else if(are(list(object@x,object@y),"xprnSet")){zx<-as.xprnSetx(object@x);zy<-as.xprnSetx(object@y)}
	 eval(parse(text=paste(class(object),'(x=zx,y=zy)',sep='')))})



setMethod("nvar", signature(x="prnSetPair"),function(x) {nrow(exprs(x@x))})
setMethod("nvar", signature(x="prnSet"),function(x) {nrow(exprs(x))})
setMethod("nsam", signature(x="prnSetPair"),function(x) {list(x=ncol(exprs(x@x)),y=ncol(exprs(x@y)))})
setMethod("nsam", signature(x="prnSet"),function(x) {ncol(exprs(x))})

#--

prep.2matrices<-function(x,y,paired=FALSE,rm.na=T){     
	if(Is(x,d1class)){x<-as.rowmatrix(x)}
	if(Is(y,d1class)){y<-as.rowmatrix(y)}
	assert.are(list(x,y),'matrix')
	nmvarx<-c('X','I')
	if(paired){nmvary<-c('X','I')} else {nmvary<-c('X','S')}
	x<-MakeNames(x,nmvar=nmvarx)                 
	if(rm.na){x<-removeNA.from.matrix(x=x,indx=FALSE)} 
	if(!is.vide(y)){
		y<-MakeNames(y,nmvar=nmvary)
		if(rm.na){y<-removeNA.from.matrix(x=y,indx=FALSE)}
		indx<-k.indSort(x=rownames(y),ref=rownames(x),ind=T,inter=T)#,rm.na=T,inter=T
		if(length(indx)==0){stop('different rownames (feature names) in x and y')}
		y<-y[indx$x,,drop=FALSE];x<-x[indx$ref,,drop=FALSE]
		if(setequal(colnames(y),colnames(x))||paired==T){
			indx<-k.indSort(x=colnames(y),ref=colnames(x),ind=T,inter=T)
			if(is.vide(indx)){stop('no paired x and y')}
			y<-y[,indx$x,drop=FALSE];x<-x[,indx$ref,drop=FALSE]
			stopifnot(identical(colnames(x),colnames(y)))
		}
		}
	stopifnot(identical(rownames(x),rownames(y)))
	list(x=x,y=y,info=list(paired=paired))
}

k.convert.paired2simple<-function(x,y){
	assert.are(list(x,y),'matrix')
	
	z<-prep.2matrices(x=x,y=y,paired=T)
	z$x-z$y
}

k.new.dataset<-function( x,rm.na=T,pdata=NULL,fdata=NULL,nmvar=c('X','I'),nmvar.ini=c(1,1)){
	
	if(Is(x,d1class)){x<-as.rowmatrix(x)}
	assert.is(x,d2class)
	assert.are(list(pdata,fdata),c(d1class,d2class,'NULL'))
	x<-MakeNames(x=x)
	
	if(length(fdata)>0){
		if(Is(fdata,d1class)){
		    if(length(fdata)==1){fdata<-rep(fdata,nrow(x))}
		    fdata<-as.colmatrix(fdata)}
		else if(Is(fdata,d2class)){if(nrow(fdata)==1){fdata<-t(fdata)}}
		
		assert.is(fdata,d2class)
		stopifnot(nrow(fdata)>=nrow(x))
		
		if(!is(fdata,'data.frame')){
			fdata<-as.data.frame(fdata,row.names=rownames(fdata),stringsAsFactors=FALSE)}
		
		num.nmp<-suppressWarnings(as.numeric(rownames(fdata)))
		if(is(fdata,'data.frame')&length(num.nmp)>0&all(is(num.nmp,'numeric'))){
			fdata<-fdata[1:nrow(x),,drop=FALSE]
			rownames(fdata)<-rownames(x)}
	}
	if(length(pdata)>0){
		if(Is(pdata,d1class)){
		    if(length(pdata)==1){pdata<-rep(pdata,ncol(x))}
		    pdata<-as.colmatrix(pdata)}
		else if(Is(pdata,d2class)){if(nrow(pdata)==1){pdata<-t(pdata)}}
		assert.is(pdata,d2class)
		stopifnot(nrow(pdata)>=ncol(x))
		
		if(!is(pdata,'data.frame')){
			pdata<-as.data.frame(pdata,row.names=rownames(pdata),stringsAsFactors=FALSE)}
		num.nmp<-suppressWarnings(as.numeric(rownames(pdata)))
		if(is(pdata,'data.frame')&length(num.nmp)>0&all(is(num.nmp,'numeric'))){
			pdata<-pdata[1:ncol(x),,drop=FALSE]
			rownames(pdata)<-colnames(x)}
	}
    
	
	if(rm.na){x<-removeNA.from.matrix(x=x,indx=FALSE)}

    if(length(fdata)>0){

        indx<-indSortAsY(x=rownames(fdata),y=rownames(x),ind=T,inter=T)$x
	fdata<-fdata[indx,,drop=FALSE]

        stopifnot(identical(rownames(fdata),rownames(x)))
	stopifnot(nrow(fdata)>=nrow(x))
	}
    else{fdata<-data.frame()}
    
    if(length(pdata)>0){
	
        indx<-indSortAsY(x=rownames(pdata),y=colnames(x),inter=T,ind=T)$x
	pdata<-pdata[indx,,drop=FALSE]
        
        stopifnot(identical(rownames(pdata),colnames(x)))
	stopifnot(nrow(pdata)>=ncol(x))

	}
    else{pdata<-data.frame()}

	#---------

	
    list(x=x,pdata=pdata,fdata=fdata,positive=is.positive(x))}




k.new.dataset.pair<-function(x,y,paired=FALSE,rm.na=T,x.pdata=NULL,y.pdata=NULL,x.fdata=NULL,y.fdata=NULL,nmvar=c('X','I','S'),nmvar.ini=c(1,1)){

	assert.are(list(x.pdata,y.pdata,x.fdata,y.fdata),c(d1class,d2class,'NULL'))
	
	stopifnot(length(nmvar)==3 & length(nmvar.ini)==2)
	if(paired){nmvarx<-nmvary<-nmvar[1:2]}
	else{nmvarx<-nmvar[1:2];nmvary<-nmvar[c(1,3)]}
	z<-prep.2matrices(x=x,y=y,paired=paired,rm.na=rm.na)
	x<-z$x;y<-z$y
	zx<-k.new.dataset(x=x,rm.na=rm.na,pdata=x.pdata,fdata=x.fdata,nmvar=nmvarx,nmvar.ini=nmvar.ini)
	zy<-k.new.dataset(x=y,rm.na=rm.na,pdata=y.pdata,fdata=y.fdata,nmvar=nmvary,nmvar.ini=nmvar.ini)

	xx<-zx$x;yy<-zy$x
	x.pdata<-zx$pdata
	y.pdata<-zy$pdata
	x.fdata<-zx$fdata
	y.fdata<-zy$fdata
               
	nids.f<-!identical(x.fdata,y.fdata)
	ids.f<-(identical(rownames(x.fdata),rownames(xx)))|(identical(rownames(y.fdata),rownames(yy)))
	if(nids.f & ids.f){
		fdata<-cbindx(x.fdata,y.fdata)
		colnames(fdata)<-make.unique(colnames(fdata))}
	else if(identical(x.fdata,y.fdata)|| is.vide(y.fdata)){fdata<-x.fdata}
	else if(is.vide(x.fdata)){fdata<-y.fdata}
	
	#----------

	pos<-is.positive(xx)&is.positive(yy)
	#-----------
	if(paired){stopifnot(identical(colnames(xx),colnames(yy)))}

	stopifnot(identical(rownames(xx),rownames(yy)))
	if(nrow(x.pdata)>0){stopifnot(identical(colnames(xx),rownames(x.pdata)[1:ncol(xx)]))}
	if(nrow(y.pdata)>0){stopifnot(identical(colnames(yy),rownames(y.pdata)[1:ncol(yy)]))}
	if(nrow(fdata)>0){stopifnot(identical(rownames(fdata)[1:nrow(xx)],rownames(xx)))}
	list(x=xx,y=yy,x.pdata=x.pdata,y.pdata=y.pdata,fdata=fdata,positive=pos,paired=paired)
}
#----
k.matrix2ExpressionSet<-function(x,rm.na=T,pdata=NULL,fdata=NULL,annot=NULL,...){nmvar<-c('X','I');nmvar.ini<-c(1,1)
	z<-k.new.dataset(x=x,rm.na=rm.na,pdata=pdata ,fdata=fdata,nmvar=nmvar,nmvar.ini=nmvar.ini)
	xx<-z$x;pdata<-z$pdata;fdata<-z$fdata
	if(length(annot)==0){annot<-''}
	if(!is.character(annot)){annot<-as.character(annot)}
	#-------
	if(all(is.vide(pdata))|nrow(pdata)==0){
		pdata <- data.frame(dummy = factor(rep(0, ncol(xx))))
		rownames(pdata) <- colnames(xx)}
	#-------

	pdata<-pdata[1:ncol(xx),,drop=FALSE]
	if(all(is.vide(fdata))|nrow(fdata)==0){
		stopifnot(identical(colnames(xx),rownames(pdata)[1:ncol(xx)]))
	ed <- new("ExpressionSet", phenoData = as(pdata, "AnnotatedDataFrame"),exprs = xx,annotation=annot,...)}
	else{fdata<-fdata[1:nrow(xx),,drop=FALSE]
		stopifnot(identical(rownames(xx),rownames(fdata)[1:nrow(xx)]))
	ed <- new("ExpressionSet", phenoData = as(pdata, "AnnotatedDataFrame"),featureData=as(fdata, "AnnotatedDataFrame"),
			exprs = xx,annotation=annot, ...)}
        ed
}

k.matrix2ExpressionSetPair<-function(x,y,paired=FALSE,rm.na=T,x.pdata=NULL,y.pdata=NULL,x.fdata=NULL,y.fdata=NULL,annot=c('case','control'),...){
	nmvar<-c('X','I','S');nmvar.ini<-c(1,1)
	
	z<-k.new.dataset.pair(x=x,y=y,paired=paired,x.pdata=x.pdata,y.pdata=y.pdata,rm.na=rm.na,
			x.fdata=x.fdata,y.fdata=y.fdata,nmvar=nmvar,nmvar.ini=nmvar.ini)

	if(is.vide(annot)){annot<-c('','','')}
	if(length(annot)==1){annot<-rep(annot,3)}
	if(length(annot)==2){annot<-c(annot,'')}
	zx<-k.matrix2ExpressionSet(x=z$x,rm.na=rm.na,pdata=z$x.pdata,fdata=z$fdata,annot=annot[1],...)
	zy<-k.matrix2ExpressionSet(x=z$y,rm.na=rm.na,pdata=z$y.pdata,fdata=NULL,annot=annot[2],...)
	list(x=zx,y=zy,annotation=annot)
}



k.xprnSet<-function(phenoData=NULL, exprs, featureData=NULL,annotation=NULL,...){
	ed<-k.matrix2ExpressionSet(x=exprs,rm.na=T,pdata=phenoData,fdata=featureData,annot=annotation,...)
	new("xprnSet", es = ed)
}




#---
k.prep.xyfData<-function(x,y,x.fdata=NULL,y.fdata=NULL,fdata=NULL){
	
	ff<-xf<-yf<-NULL
	if(is(x,'prnSet')){
		xf<-try(fData(x),silent=T);x<-exprs(x);
		if(!is.err(xf)){xf<-k.prep.fData(exprs=x,featureData=xf)}
			else{xf<-NULL}
		}
	if(is(y,'prnSet')){
		yf<-try(fData(y),silent=T);y<-exprs(y)
		if(!is.err(yf)){yf<-k.prep.fData(exprs=y,featureData=yf)}
			else{yf<-NULL}
		}
	assert.are(list(x,y),c('matrix','array'))
	x<-MakeNames(x);y<-MakeNames(y)
	fint<-intersect(rownames(x),rownames(y))
	
	if(!is.vide(fdata)){ff<-k.prep.fData(exprs=x[fint,,drop=F],featureData=fdata)}
		else{ff<-NULL}
	
	nxf<-k.prep.fData(exprs=x,featureData=x.fdata)
	nyf<-k.prep.fData(exprs=y,featureData=y.fdata)
	
	list(x.fdata=cbindx(ff,cbindx(xf,nxf)),y.fdata=cbindx(yf,nyf))

}
k.prep.fData<-function(exprs,featureData=NULL){
	if(is.Matrix(exprs)){exprs<-as.matrices(exprs)}
		
	exprs<-MakeNames(exprs);assert.is(exprs,'matrix')
	if(length(featureData)==0){featureData<-NULL}
		
	if(Is(featureData,d1class)){
		if(length(featureData)==1){
			featureData<-rep(featureData,nrow(exprs))
			names(featureData)<-rownames(exprs)
		}
		featureData<-as.colmatrix(featureData)
		}
	if(Is(featureData,'matrix')){
		featureData<-as.data.frame(featureData,stringsAsFactors =F)
		#featureData<-as.data.frame(featureData,row.names=rownames(exprs))
		featureData<-MakeColNames(featureData,nmvar='FALSE')}
	
	assert.is(featureData,c('data.frame','NULL'))
	
	if(Is(featureData,'data.frame')){
		num.nmp<-suppressWarnings(as.numeric(rownames(featureData)))
		if(length(num.nmp)>0&all(is.finite(num.nmp))){
			featureData<-featureData[1:nrow(exprs),,drop=FALSE]
			rownames(featureData)<-rownames(exprs)}
		featureData<-featureData[rownames(exprs) ,,drop=F]
		stopifnot(nrow(exprs) == nrow(featureData) && all(rownames(exprs) == rownames(featureData)))
		}
	
	assert.is(featureData,c('data.frame','NULL'))
	featureData	
}
k.prep.featureData<-function(exprs,featureData){
	if(Is(featureData,'missing')){featureData<-NULL}
	featureData<-k.prep.fData(exprs=exprs,featureData=featureData)
	
	if(Is(featureData,'data.frame')){featureData <- as(featureData, "AnnotatedDataFrame")}
	assert.is(featureData,c('AnnotatedDataFrame','NULL'))
	featureData	
}
k.prep.phenoData<-function(exprs,phenoData){
	if(is.Matrix(exprs)){exprs<-as.matrices(exprs)}
		
	exprs<-MakeNames(exprs);assert.is(exprs,'matrix')
	if(Is(phenoData,'missing')){phenoData<-NULL}
	if(length(phenoData)==0){phenoData<-NULL}
		
	if(Is(phenoData,d1class)){
		if(length(phenoData)==1){
			phenoData<-rep(phenoData,ncol(exprs))
			names(phenoData)<-colnames(exprs)
		}
		phenoData<-as.colmatrix(phenoData)
		}
	if(Is(phenoData,'matrix')){
		phenoData<-as.data.frame(phenoData,stringsAsFactors =F)
		#phenoData<-data.framex(phenoData,row.names=colnames(exprs))
		phenoData<-MakeColNames(phenoData,nmvar='FALSE')}
	
	assert.is(phenoData,c('data.frame','NULL',"AnnotatedDataFrame"))
	
	if(is.vide(phenoData)){
		phenoData <- data.frame(dummy = factor(rep(0, ncol(exprs))))
		rownames(phenoData) <- colnames(exprs)
	}
	
	if(Is(phenoData,'data.frame')){
		num.nmp<-suppressWarnings(as.numeric(rownames(phenoData)))
		if(length(num.nmp)>0&all(is.finite(num.nmp))){
			phenoData<-phenoData[1:ncol(exprs),,drop=FALSE]
			rownames(phenoData)<-colnames(exprs)}
		phenoData<-phenoData[colnames(exprs),,drop=F]
		stopifnot(ncol(exprs) == nrow(phenoData) && all(colnames(exprs) == rownames(phenoData)))
		phenoData <- as(phenoData, "AnnotatedDataFrame")}
	
	
	assert.is(phenoData,c('AnnotatedDataFrame'))
	phenoData	
}

#---
setMethod("xprnSet", signature(phenoData = "all.matrixORnumeric", exprs = "missingOrNULL",featureData = "missingOrNULL"), function(phenoData, exprs,featureData, ...)
{xprnSet(exprs = phenoData, ...)})
setMethod("XprnSet", signature(phenoData = "all.matrixORnumeric", exprs = "missingOrNULL",featureData = "missingOrNULL"), function(phenoData, exprs,featureData, ...)
{	exprs<-phenoData;exprs<-MakeNames(exprs)
	if(!all(as.numeric(exprs) > 0, na.rm = TRUE))
	{ message("XprnSet requires all positive; consider xprnSet")}
	if(nis.Matrix(exprs)){exprs<-Matrix(exprs)}
	new("XprnSet", xprnSet(exprs = exprs, ...))
})




setMethod("xprnSet", signature(phenoData = "all.pfdata", exprs = "all.matrixORnumeric", featureData="all.pfdata"), function(phenoData, exprs,featureData, ...){
	exprs<-MakeNames(exprs)
	if(class(exprs)%in%"Numeric"){exprs<-Matrix(as.rowmatrix(exprs))}
	if(is(exprs,"numeric")){exprs<-as.rowmatrix(exprs)}

	if(missing(phenoData)){phenoData<-NULL}
	if(missing(featureData)){featureData<-NULL}
	featureData<-k.prep.featureData(exprs=exprs,featureData=featureData)
	phenoData<-k.prep.phenoData(exprs=exprs,phenoData=phenoData)
	
	if(is.vide(featureData)){xprnSet(phenoData = phenoData, exprs = exprs,...)}
	else{xprnSet(phenoData = phenoData, exprs = exprs,featureData = featureData, ...)}
	#k.xprnSet(phenoData = phenoData, exprs = exprs,featureData=featureData,annot=annotation,...)
})
setMethod("XprnSet", signature(phenoData = "all.pfdata", exprs = "all.matrixORnumeric", featureData="all.pfdata"), function(phenoData, exprs,featureData, ...){
	if(!all(as.numeric(exprs) > 0, na.rm = TRUE)){ message("XprnSet requires all positive; consider xprnSet")}
	exprs<-MakeNames(exprs)
	if(is(exprs,"numeric")){exprs<-as.rowmatrix(exprs)}
	if(nis.Matrix(exprs)){exprs<-Matrix(exprs)}

	if(missing(phenoData)){phenoData<-NULL}
	if(missing(featureData)){featureData<-NULL}
	featureData<-k.prep.featureData(exprs=exprs,featureData=featureData)
	phenoData<-k.prep.phenoData(exprs=exprs,phenoData=phenoData)
	#k.xprnSet(phenoData = phenoData, exprs = exprs,featureData=featureData,annot=annotation,...)
	zo<-xprnSet(exprs = exprs,phenoData = phenoData,featureData = featureData,...)
	new("XprnSet",zo)
})




#
setMethod("xprnSubset", signature(object = "XprnSet"), function(object, ...){
	object@es <- xprnSubset(object = object@es, ...)
	object
})
setMethod("xprnSetPair", signature(x = "XprnSet", y = "XprnSet", factor.name = "missing"), function(x, y, factor.name)
{
	xprnSetPair(x = as.xprnSetx(x,base=exp(1)), y = as.xprnSetx(y,base=exp(1)))
})
##--------------------------------------
setMethod("xprnSetPair", signature(x = "all.matrixORnumeric", y = "all.matrixORnumeric", factor.name = "missing"), function(x, y, factor.name)
{	x<-MakeNames(x);y<-MakeNames(y)
	if(is.Matrix(x)){x<-as.matrices(x)}
	if(is.Matrix(y)){y<-as.matrices(y)}
		x<-xprnSet(exprs=x)
		y<-xprnSet(exprs=y)
	new("xprnSetPair", x = x, y = y)
		
})
#
setMethod("logb", signature(x = "XprnSet"), function(x,base=exp(1)){
	as.xprnSetx(x,base=base)
})
setMethod("exp", signature(x = "xprnSet"), function(x){ base<-exp(1)
	as.XprnSetx(x,base=base)
})
setMethod("log", signature(x = "ExpressionSet"), function(x,base=exp(1)){
	if(base!=exp(1)){z<-log(exprs(x),base=base)}
	else{z<-log(exprs(x))}
	exprs(x)<-z;x
})
setMethod("logb", signature(x = "ExpressionSet"), function(x,base=exp(1)){
	z<-logb(exprs(x),base=base)
	exprs(x)<-z;x
})
setMethod("log2", signature(x = "ExpressionSet"), function(x){
	z<-log2(exprs(x))
	exprs(x)<-z;x
})

#-----

setAs(from = "XprnSet", to = "ExpressionSet", function(from){from@es})
setAs(from = "xprnSet", to = "ExpressionSet", function(from){from@es})
setAs(from = "ExpressionSet", to = "XprnSet", function(from){new("XprnSet", es = from)})
setAs(from = "ExpressionSet", to = "xprnSet", function(from){new("xprnSet", es = from)})
#

nXprnSet<-function(x, pdata = NULL, fdata = NULL, annot =
                 NULL,...){z<-k.matrix2ExpressionSet(x=x,rm.na=TRUE,pdata=pdata,fdata=fdata,annot=annot,...); as.XprnSetx(object=z)}
nxprnSet<-function(x, pdata = NULL, fdata = NULL, annot =
                 NULL,...){z<-k.matrix2ExpressionSet(x=x,rm.na=TRUE,pdata=pdata,fdata=fdata,annot=annot,...); as.xprnSetx(object=z)}
nxprnSetPair<-function(x,y=NULL,paired=FALSE,x.pdata=NULL,y.pdata=NULL,x.fdata=NULL,y.fdata=NULL,fdata=NULL,annot=c('case','control'),factor.name,...){rm.na<-TRUE

#if(!is.vide(fdata)){x.fdata<-cbind(fdata,x.fdata,y.fdata);y.fdata<-NULL}
	if(are(list(x,y),c('numeric','matrix'))){
		if(is.Matrix(x)){x<-as.matrices(x)}
		if(is.Matrix(y)){y<-as.matrices(y)}
		fx<-k.prep.xyfData(x=x,y=y,x.fdata=x.fdata,y.fdata=y.fdata,fdata=fdata)	
		z<-k.matrix2ExpressionSetPair(x=x,y=y,paired=paired,rm.na=rm.na,x.pdata=x.pdata,y.pdata=y.pdata,x.fdata=fx$x.fdata,y.fdata=fx$y.fdata,annot=annot,...)
		x<-as.xprnSetx(object=z$x,base=NULL);y<-as.xprnSetx(object=z$y,base=NULL)
		}
	if(are(list(x,y),'xprnSet')){
		fx<-k.prep.xyfData(x=x,y=y,x.fdata=x.fdata,y.fdata=y.fdata,fdata=fdata)		
		z<-k.matrix2ExpressionSetPair(x=exprs(x),y=exprs(y),paired=FALSE,rm.na=rm.na,
			x.pdata=pData(x),y.pdata=pData(y),x.fdata=fx$x.fdata,y.fdata=fx$y.fdata,
			annot=paste(annotation(x),' vs. ',annotation(y),sep=''),...)
		
		if(are(list(x,y),'XprnSet')){message('making log of x and y (class XprnSet)')
			xx<-as.xprnSetx(object=log(z$x),base=NULL)
			yy<-as.xprnSetx(object=log(z$y),base=NULL)
			}
		else{
			xx<-as.xprnSetx(object=z$x,base=NULL)
			yy<-as.xprnSetx(object=z$y,base=NULL)	
		}
		
		zo<-xprnSetPair(x=xx,y=yy)
	}
	else if(is(x,"XprnSet")&&(missing(y)||is.vide(y))){
		x<-as.xprnSetx(object=x,base=exp(1))
		zo<-xprnSetPair(x=x,factor.name=factor.name)
	}
	else if(is(x,"xprnSet")&&(missing(y)||is.vide(y))){
		zo<-xprnSetPair(x=x,factor.name=factor.name)
	}
	else{stop('bad inputs')}
	zo
	
}



#--------------
x2X<-function(x,base=NULL){
	if(is(x,'matrix')){zm<-x}
	else if(is(x,"ExpressionSet")){zm<-exprs(x)}
	else{stop('bad input')}
	
	if(is.vide(base)){xm<-zm}
	else{
		if(base!=exp(1)){xm<-base^(zm)}
		else{xm<-exp(zm)}
		}
	
	if(nis.Matrix(x)){xm<-Matrix(xm)}
	if(Is(x,c('matrix','numeric'))){zo<-xm}
	else if(is(x,"ExpressionSet")){
		exprs(x)<-xm
		zo<-x}
	
zo	
}
X2x<-function(x,base=NULL){
	if(is.Matrix(x)){zm<-as.matrices(x)}
	else if(nis.Matrix(x)){zm<-x}
	else if(is(x,"ExpressionSet")){zm<-exprs(x)}
	else{stop('bad input')}
	
	if(is.vide(base)){xm<-zm}
	else{
		if(base!=exp(1)){xm<-log(x=zm, base = base)}
		else{xm<-log(x=zm);base<-'e'}
	}
	
	if(Is(x,c('matrix','numeric'))){zo<-xm}
	else if(is(x,"ExpressionSet")){
		exprs(x)<-xm
		zo<-x}
zo	
}
as.xprnSetx <- function(object,base=exp(1)){
	if(is.xprnSet(object)){zo<-object}
	else if(nis.Matrix(object)){zo<-xprnSet(exprs = object)}
	else if(is.Matrix(object)){
		zm<-X2x(x=object,base=base)
		if(is.Matrix(zm)){zm<-as.matrices(zm)}
		zo<-xprnSet(exprs = zm)}
	else if(is.ExpressionSet(object)){zo<-new("xprnSet", es = object)}
	else if(is.XprnSet(object)){
		zm<-X2x(x=exprs(object),base=base)
                if(is.Matrix(zm)){zm<-as.matrices(zm)}
		z<-as(object, "ExpressionSet")
		exprs(z)<-zm;zo<-new("xprnSet", es = z)}
	else if(is.prnSetPair(object)){
		z<-k.convert.paired2simple(x=exprs(object@x),y=exprs(object@y))
		if(is.XprnSetPair(object)){
			z<-X2x(x=z,base=base)
			if(is.Matrix(z)){z<-as.matrices(z)}}
		zo<-xprnSet(phenoData = cbindx(pData(object@x),pData(object@y)), exprs = z, featureData=fData(object))
	}
		
	validObject(zo);zo	
}

as.XprnSetx <- function(object,base=exp(1)){
	if(is.XprnSet(object)){zo<-object}
	else if(is.Matrix(object)){zo<-XprnSet(exprs = object)}
	else if(nis.Matrix(object)){
		zm<-x2X(x=object,base=base)
		if(!is.Matrix(zm)){zm<-Matrix(zm)}
		zo<-XprnSet(exprs = zm)}
	else if(is.ExpressionSet(object)){zo<-new("XprnSet", es = object)}
	else if(is.xprnSet(object)){
		zm<-x2X(x=exprs(object),base=base)
                if(nis.Matrix(zm)){zm<-Matrix(zm)}
		z<-as(object, "ExpressionSet")
		exprs(z)<-zm;zo<-new("XprnSet", es = z)}
	else if(is.prnSetPair(object)){
		z<-k.convert.paired2simple(x=exprs(object@x),y=exprs(object@y))
		if(is.xprnSetPair(object)){
			z<-x2X(x=z,base=base)
			if(nis.Matrix(z)){z<-Matrix(z)}}
		zo<-XprnSet(phenoData = cbindx(pData(object@x),pData(object@y)), exprs = z, featureData=fData(object))
	}
		
	validObject(zo);zo

}




as.xprnSetPairx<- function(object,base=NULL){
	if(is.xprnSetPair(object)){zo<-object}
	else if(is.XprnSetPair(object)){
		xx<-X2x(exprs(object)$x,base=base)
		yy<-X2x(exprs(object)$y,base=base)
		zo<-nxprnSetPair(x=xx,y=yy,x.pdata=pData(object@x),y.pdata=pData(object@y),fdata=fData(object),annot=c(annotation(object@x),annotation(object@x)))}
	
	validObject(zo);zo
}



#
prnSet2matrix<-function(x,y=NULL,paired=FALSE){
	assert.is(x,c('xprnSet','xprnSetPair','matrix'))
	assert.is(y,c('xprnSet','NULL','matrix'))
	if(is(x,'prnSet')){x<-exprs(x)}
	if(is(y,'prnSet')){y<-exprs(y)}
	if(is(x,'prnSetPair')){y<-exprs(x)$y;x<-exprs(x)$x}
	y<-MakeNames(y);x<-MakeNames(x)

	if(!is.vide(y)&&paired==T){
		#z<-prep.2matrices(x=x,y=y,paired=paired,rm.na=rm.na)
		#x<-z$x;y<-z$y;paired<-z$info$paired
		pok<-is.paired(x=x,y=y)
		if(!pok&&paired==T){stop('data is not paired')}
		if(paired==T){x<-x-y;y<-NULL;paired<-FALSE}}
	assert.is(x,'matrix')
	assert.is(y,c('matrix','NULL'))
	
	list(x=x,y=y,paired=paired)}

#============pvalues=======================

k.get.from.funtest<-function(x,y=NULL,fun=t.test,opt='p.value',...){rm.na<-FALSE
	x<-as.numeric(x);if(!is.vide(y)){y<-as.numeric(y)}
	assert.are(list(x,y),c('numeric','NULL'))

	#if(rm.na){x<-x[!is.na(x)];y<-y[!is.na(x)]}
		x<-x[!is.na(x)]
	if(is.vide(y)){eval(parse(text=paste('zo<-try(fun(x=x,...)$',opt,')',sep='')))}
	else{y<-y[!is.na(y)];eval(parse(text=paste('zo<-try(fun(x=x,y=y,...)$',opt,')',sep='')))}
	if(is.err(zo)){zo<-as.numeric(NA)}
			zo
}

k.get.from.funtest.mat<-function(x,y=NULL,fun=t.test,paired=FALSE,opt='p.value',...){rm.na<-FALSE
	
	z<-prnSet2matrix(x=x,y=y,paired=paired);x<-z$x;y<-z$y
	assert.is(x,'matrix');assert.is(y,c('matrix','NULL'))

	argk<-uniquex(list(...))
	
	zo<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){
		xx<-x[i,];yy<-y[i,]
		if(paired && identical(names(xx),names(yy))){
			ix<-intersect(which(!is.na(xx)),which(!is.na(xx)))
			if(length(ix)==0){warning('non finite elements');return(as.numeric(NA))}
			xx<-xx[ix];yy<-yy[ix]
			}
		
		do.call(k.get.from.funtest,c(list(x=xx,y=yy,fun=fun,rm.na=rm.na,opt=opt),argk))[1]
		})
	names(zo)<-rownames(x)
zo
}

get.other.from.testfun<-function(x,y=NULL,fun=t.test,paired=FALSE,opt='parameter',...){
	k.get.from.funtest.mat(opt=opt,x=x,y=y,fun=fun,paired=paired,...)
}
get.pvalues<-function(x,y=NULL,fun=t.test,paired=FALSE,...){
	zo<-k.get.from.funtest.mat(opt='p.value',x=x,y=y,fun=fun,paired=paired,...)
	stopifnot(all(is.prob(zo[is.finite(zo)])))
	zo
}
get.stats<-function(x,y=NULL,fun=t.test,paired=FALSE,...){
	k.get.from.funtest.mat(opt='statistic',x=x,y=y,fun=fun,paired=paired,...)
}
get.df<-function(x,y=NULL,paired=FALSE,var.equal=T,fun=t.test,...){
	if(is(x,'numeric')){x<-as.rowmatrix(x)}
	if(is(y,'numeric')){y<-as.rowmatrix(y)}
	dfx<-get.other.from.testfun(x=x,y=y,fun=fun,paired=paired,opt='parameter',var.equal=var.equal,...)

dfx
}
stat2pval<-function(x,pFUN=pt,sym.distrib=T,alternative='two.sided',...){#alternative = "greater" is the alternative that x has a larger mean than y
	assert.is(x, "numeric")
	pnm<-getName(pFUN)
	
	if(sym.distrib&&pnm%in%c('pchisq')&&alternative%in%c('t','two.sided')){warning('setting sym.distrib=FALSE');sym.distrib<-FALSE}
	if(!sym.distrib&&pnm%in%c('pt')&&alternative%in%c('t','two.sided')){warning('setting sym.distrib=T');sym.distrib<-T}
	
	if(alternative%in%c('g','greater','t','two.sided')){lower.tail<-FALSE}
	else if(alternative%in%c('l','less')){lower.tail<-T}
	else{stop('bad alternative')}
	pv<-pFUN(q=x,lower.tail=lower.tail,...)
	if(alternative%in%c('t','two.sided')&&sym.distrib){pv<-2*pmin(pv,pFUN(q=x,lower.tail=!lower.tail,...))}
	else if(alternative%in%c('t','two.sided')&&!sym.distrib){
		pv<-pFUN(q=x,lower.tail=FALSE,...)+pFUN(q=-x,lower.tail=T,...)}
		pv<-pmin(pv,1)

	stopifnot(all(is.prob(pv[is.finite(pv)])))
	names(pv)<-names(x);pv
}
pval2stat<-function(x,qFUN=qt,alternative='two.sided',...){
	assert.is(x, "numeric")
	stopifnot(all(is.prob(x[is.finite(x)])))
	if(alternative%in%c('g','greater','t','two.sided')){lower.tail<-FALSE}
	else if(alternative%in%c('l','less')){lower.tail<-T}
	else{stop('bad alternative')}
	if(alternative%in%c('t','two.sided')){warning('wrong statistics from two.sided p-values');x<-x/2}
	
	st<-qFUN(p=x,lower.tail=lower.tail,...)
	names(st)<-names(x);st
}
#
data2statANDpvalue<-function(x,y=NULL,fun=abs.t.test,paired=FALSE,opt=c('both','all','b','a','pvalue','p','stat','s'),...){opt<-opt[1]
  assert.is(x,c('prnSet','prnSetPair','matrix','numeric'))
  assert.is(y,c('prnSet','NULL','matrix','numeric'))
  spfun<-As.function(fun) 
  stopifnot(is(spfun,'function'))
  
  arglis.fun<-uniquex(c(list(paired=paired),list(...)))
		
        infox<-list()
    if(print.info){
      info.fun<-get.args(fun=fun,...)
      info.fun<-uniquex(c(arglis.fun,info.fun))
      infox<-info.fun[!names(info.fun)%in%c('x','y')]
    }
    stat<-pval<-NULL
if(opt%in%c('stat','s','both','all','b','a')){
	stat<-do.call(get.stats,c(list(x=x,y=y,fun=fun),arglis.fun))}
if(opt%in%c('pvalue','p','both','all','b','a')){
	pval<-do.call(get.pvalues,c(list(x=x,y=y,fun=fun),arglis.fun))}

  list(stat=stat,pvalue=pval,info=infox)
}


stat.toANDfrom.pvalue<-function(stat=NULL,pvalue=NULL,cFUN,alternative="two.sided",sym.distrib=T,...){
	stopifnot(any.mustbe.nonnull(pvalue,stat))
	assert.is(cFUN,c('character','function'))
	if(is.vide(pvalue)){
		assert.is(stat,'numeric')
		stat<-MakeNames(stat)
		stopifnot(any(is.finite(stat)))
		orig.stat<-stat}
	if(is.vide(stat)){
		assert.is(pvalue,'numeric')
		pvalue<-MakeNames(pvalue)
		stopifnot(any(is.finite(pvalue)))
		stopifnot(all(is.prob(pvalue[is.finite(pvalue)])))
		orig.pvalue<-pvalue}
		
	if(length(pvalue)==0){
		pvalue<-do.call(stat2pval,uniquex(c(list(x=stat,pFUN=cFUN,alternative=alternative,sym.distrib=sym.distrib),list(...))))}
	if(length(stat)==0){
		stat<-do.call(pval2stat,uniquex(c(list(x=pvalue,qFUN=cFUN,alternative=alternative),list(...))))}
		
	stopifnot(length(pvalue)==length(stat) && identical(names(pvalue),names(stat)))
	stopifnot(all(is.prob(pvalue[is.finite(pvalue)])))
	list(pvalue=pvalue,stat=stat,info=list(cFUN=getName(cFUN),arglis.cFUN=list(...),alternative=alternative))	
}
#===========================
nplot<-function(x,y,...){}
   
#-----------------------------
#-------------
q.absTd<-function (p, df, ncp = 0, lower.tail=FALSE,...) {
	stopifnot(all(is.prob(p[is.finite(p)])))
    at <- absTd(df = df, ncp = ncp)
    q.r(at)(p = p,lower.tail=lower.tail, ...)
};attr(q.absTd,'name')<-'q.absTd'
d.absTd<-function (x, df, ncp = 0, ...) {
	stopifnot(all(x>=0))
    at <- absTd(df = df, ncp = ncp)
    d(at)(x = x, ...)
};attr(d.absTd,'name')<-'d.absTd'
p.absTd<-function (q, df, ncp = 0, lower.tail=FALSE,...) {
	ix<-which(q<0)
	q<-abs(q)
	#stopifnot(all(q[is.finite(q)]>=0))
    at <- absTd(df = df, ncp = ncp)
    px<-p(at)(q = q, lower.tail=lower.tail,...)
    if(lower.tail&&length(ix)>0){px[ix]<-0}
else if(!lower.tail&&length(ix)>0){px[ix]<-1}
px
};attr(p.absTd,'name')<-'p.absTd'
dabsTd<-function(x,df,ncp=0){dt(x=x,df=df,ncp=ncp)+dt(x=-x,df=df,ncp=ncp)}
attr(dabsTd,'name')<-'dabsTd'
pabsTd<-function(q,df,ncp=0,lower.tail =T){p<-2*pt(q=q,df=df,ncp=ncp,lower.tail = FALSE);if(lower.tail){p<-1-p};p}
attr(pabsTd,'name')<-'pabsTd'
qabsTd<-function(p,df,ncp=0,lower.tail =T){stopifnot(all(is.prob(p) ));q<-qt(p=p/2,df=df,ncp=ncp,lower.tail = FALSE);if(lower.tail){q<--q};abs(q)}
attr(qabsTd,'name')<-'qabsTd'
rabsTd<-function(...){z<-rt(...);abs(z)};attr(rabsTd,'name')<-'rabsTd'
abs.t.test<-function(...){z<-t.test(...);z$statistic<-abs(z$statistic);z};attr(abs.t.test,'name')<-'abs.t.test'

qabsTd<-function(p,df,ncp=0,lower.tail =T){stopifnot(all(is.prob(p) ))
if(lower.tail){np<-(1-p)/2} else {np<-p/2}
q<-qt(p=np,df=df,ncp=ncp,lower.tail = FALSE);
stopifnot(all(q[is.finite(q)]>=0));q}
attr(qabsTd,'name')<-'qabsTd'
#============for LFDRs==========
setClass("est.lfdr.pvalue", representation(LFDR.hat ="numeric",p0.hat = "numeric",pvalue="numeric",info="list"))
setValidity("est.lfdr.pvalue", function(object) {
  dx<-object@LFDR.hat
  p0<-object@p0.hat
  stx<-object@pvalue
  ix<-object@info
  nvar<-length(dx)
stopifnot(!are.null(names(dx)))
stopifnot(length(p0)==1)
stopifnot(length(stx)%in%c(0,nvar))
if(length(stx)==nvar){stopifnot(identical(names(stx),names(dx)))}
if(!are.prob(list(dx[is.finite(dx)],p0[is.finite(p0)],stx[is.finite(stx)]))){
	stop('error:LFDR.hat, pvalue or p0.hat out of [0,1]')}

})


new.est.lfdr.pvalue <- function(LFDR.hat,p0.hat,pvalue,method=NULL,info=list()){err<-FALSE
  if(!is(info,'list')){info<-list()}
  assert.is(LFDR.hat,c('NULL',"numeric","logical",'try-error'))
  assert.is(pvalue,c("numeric","logical"))
  stopifnot(length(pvalue)>=length(LFDR.hat))
  if(is.error(LFDR.hat)){err<-TRUE
    LFDR.hat<-rep(as.numeric(NA),length(pvalue));names(LFDR.hat)<-names(pvalue);method<-pastex('FAILED-',method)
    p0.hat<-as.numeric(NA)}
  nvar<-length(LFDR.hat)
  names(p0.hat)<-NULL

  
  #assert.are(list(LFDR.hat,p0.hat,ncp.hat,pvalue),c("numeric","logical"))
  #if(is(LFDR.hat,'logical')){xx<-as.numeric(LFDR.hat);names(xx)<-names(LFDR.hat);LFDR.hat<-xx}
  #if(is(p0.hat,'logical')){xx<-as.numeric(p0.hat);names(xx)<-names(p0.hat);p0.hat<-xx}
  #if(is(ncp.hat,'logical')){xx<-as.numeric(ncp.hat);names(xx)<-names(ncp.hat);ncp.hat<-xx}
  if(is(pvalue,'logical')){xx<-as.numeric(pvalue);names(xx)<-names(pvalue);pvalue<-xx}
  assert.are(list(LFDR.hat,p0.hat,pvalue),'numeric')
  stopifnot(length(p0.hat)%in%c(1,nvar))
    
  #pvalue<-MakeNames(pvalue,nmvar='X')
  #LFDR.hat<-MakeNames(LFDR.hat,nmvar='X')


  LFDR.hat<-sameAsY(x=LFDR.hat,y=pvalue)

  stopifnot(identical(names(pvalue),names(LFDR.hat)))
  nvar<-length(LFDR.hat)

stopifnot(length(pvalue)%in%c(nvar))
if(!err){p0.hat<-pmin(p0.hat,1)
	stopifnot(all(is.prob(LFDR.hat[is.finite(LFDR.hat)])) && 
	  all(is.prob(p0.hat[is.finite(p0.hat)])))
}

	stopifnot(all(is.prob(pvalue[is.finite(pvalue)])))

if(is.null(method)){method<-'lfdr.PsiHat'}
  info<-info[!names(info)%in%c('x','y','W','stat','p.value','pvalue')]
  if(!print.info){info<-list()}
  infox<-uniquex(c(list(method=toupper(method)),info))
new("est.lfdr.pvalue",LFDR.hat=LFDR.hat, p0.hat=p0.hat,pvalue=pvalue,info=infox)


}



#
##-----
setMethod("[", signature(x = "est.lfdr.pvalue", i = "ANY", j = "missing"), function(x, i, j, drop){
	z<-x
    nvar<-length(x@LFDR.hat)
	z@LFDR.hat <- x@LFDR.hat[i];names(z@LFDR.hat)<-names(x@LFDR.hat)[i]
    if(length(x@pvalue)==nvar){z@pvalue <- x@pvalue [i];names(z@pvalue)<-names(x@LFDR.hat)[i]}
	if(length(x@p0.hat)==nvar){z@p0.hat <- x@p0.hat[i];names(z@p0.hat)<-names(x@LFDR.hat)[i]}
	stopifnot(validObject(z))
	z})

setReplaceMethod("names", signature(x="est.lfdr.pvalue",value='character'),function(x, value){
  nvar<-length(x@LFDR.hat)
  stopifnot(length(value)==nvar)
  names(x@LFDR.hat)<-value
  if(length(x@pvalue)==nvar){names(x@pvalue) <- value}
  if(length(x@p0.hat)==nvar){names(x@p0.hat) <- value}
  stopifnot(validObject(x))
	x})
                 
setMethod("names", signature(x = "est.lfdr.pvalue"), function(x){names(x@LFDR.hat)})
setMethod("as.numeric", signature(x="est.lfdr.pvalue"),function(x) {x@LFDR.hat})
setMethod("length", signature(x = "est.lfdr.pvalue"), function(x){length(as.numeric(x))})
setMethod('is.prob',signature(P = "est.lfdr.pvalue" ), function(P,...){
         is.prob(P@LFDR.hat,...) && is.prob(P@p0.hat,...)&& is.prob(P@pvalue,...) })

#-------------------------------
plots.est.lfdr.pvalue<-function(x,alternative='two.sided',add=FALSE,thres=0.2,...){
  assert.is(x,c("est.lfdr.pvalue", "est.lfdr.pvalue.list"))
  if(is(x,"est.lfdr.pvalue")){x<-est2list(x)}
  pv<-x$pvalue
  z<-uniquex(c(list(...),list(sub=paste('p0 = ',round(100*x$p0.hat[1],2),'%',sep=''),xlab=paste('p-values (',alternative,')',sep=''),ylab='LFDR',main=paste('LFDR by ',get.info(x)$method,sep=''))))
  if(!add){do.call(plot,c(list(x=pv,y=get.lfdr(x)),z))}
  else{do.call(points,c(list(x=pv,y=get.lfdr(x)),z))}
abline(a=0,b=1,lty=2,col='grey')
abline(h=thres,lty=3,col='grey')
}
setMethod("nplot", signature(x = 'est.lfdr.pvalue',y='missing'), function(x,y,alternative='two.sided',add=FALSE,thres=0.2,...){
  plots.est.lfdr.pvalue(x=x,alternative=alternative,add=add,thres=thres,...)
  })  

setMethod("nplot", signature(x = 'list',y='missing'), function(x,y,alternative='two.sided',add=FALSE,thres=0.2,...){
	x<-list2est(x)
	nplot(x=x,alternative=alternative,add=add,thres=thres,...)
	
})


#as.paired<-function(x,y=NULL){
#	#if(are(list(x,y),"matrixORnumeric")){zo<-prep.2matrices(x=x,y=y,paired=T,rm.na=T)}
#	if (is.matrixORnumeric(x)&&is.matrixORnumeric(y)){
#		zo<-nxprnSetPair(x=x,y=y,paired=TRUE)}
#	else if (is.Matrix(x)&&is.Matrix(y)){
#		zo<-nxprnSetPair(x=x,y=y,paired=TRUE)}
#	else if(are(list(x,y),"prnSet")){
#		zo<-nxprnSetPair(x=x,y=y,x.fdata=fData(x),y.fdata=fData(y),x.pdata=pData(x),y.pdata=pData(y),
#				 annot=c(annotation(x),annotation(y)),paired=TRUE)}
#	else if (is(x,"xprnSetPair")){
#		zo<-nxprnSetPair(x=exprs(x)$x,y=exprs(x)$y,paired=TRUE,fdata=fData(x),x.pdata=pData(x)$x,y.pdata=pData(x)$y,
#				 annot=c(annotation(x)))}
#	else{stop('bad input types')}
#	
#}
#.onUnload <- function(libpath) {    
#    library.dynam.unload("Statomica", libpath)
#}
#.onLoad <- function(libname, pkgname)  {
#	packageStartupMessage("Loading Statomica")
#
#}
