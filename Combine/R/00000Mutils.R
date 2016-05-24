#-------------classes 1d,2d------
setClassUnion("nd1class",c("numeric","character","logical","integer"))
setClassUnion("nd2class",c("matrix","data.frame"))

setGeneric("is_positive", function(x) standardGeneric("is_positive"))
setMethod("is_positive", signature(x = "nd2class"), function(x){x<-as.numeric(x);is_positive(x)})
setMethod("is_positive", signature(x = "nd1class"), function(x){all(x[is.finite(x)]>0)})

is_nd1class<-function(x){is_any(x,c("numeric","character","logical","integer"))}
is_nd2class<-function(x){is_any(x,c("matrix","data.frame"))}
is_error <- function(object){is(object, "try-error")||is_nothing(object)}
is_unk <- function(object){all(is_any(object,c("try-error","NULL")))||all(is.na(object))||is_nothing(object)}
are_unk <- function(object){all(sapply(object, is_unk))}
is_unique<-function(object){
	if(length(object)==0) {return(T)}
	assert.is(object,"nd1class")
	length(unique(object))==length(object)
}
are_null <- function(object){are(object, "NULL")}
#------------
as_rowmatrix<-function(x){
	if (is_nd1class(x)){zo<-matrix(x,1,length(x));colnames(zo)<-names(x);rownames(zo)<-'X1'}
	else if (is_nd2class(x) && nrow(x)==1){zo<-x}
	else{zo<-x}
	zo
}
as_colmatrix<-function(x){
	if (is_nd1class(x)){zo<-matrix(x,length(x),1);rownames(zo)<-names(x);colnames(zo)<-'I1'}
	else if (is_nd2class(x) && ncol(x)==1){zo<-x}
	else{zo<-x}
	zo
}
#
k.add.NA.col.dataframe<-function(x,ind){
	if(is_nd1class(x)){zo<-x[ind];return(zo)}
	
	assert.is(x,"nd2class")
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
k.indSort<-function(x=NULL,ref=NULL,inter=F,ind=T){#work with x and ref having each non repeated elements
	if(length(ref)==0){stop("ref is empty")}
	if(!is_unique(ref) || !is_unique(x)){stop("x or ref has repeated elements")}
	assert.are(list(x,ref),"nd1class")
	
	ox<-x;oref<-ref
	names(ox)<-ox;names(oref)<-oref

	xx<-1:length(x);names(xx)<-x
	yy<-1:length(ref);names(yy)<-ref
	fint<-as.character(intersect(ref,x))
	if(inter){
		if(length(fint)==0){stop("x and ref has nothing in common")}
		xx<-xx[fint]
		yy<-yy[fint]
		if(!identical(names(xx),names(yy)) || !identical(names(xx),fint)){stop("error")}
		zz<-nms<-fint
	}
	else{
		ny.nx<-as.character(setdiff(ref,x))
		nx.ny<-as.character(setdiff(x,ref))
		nms<-zz<-c(ref,nx.ny)
		xx<-xx[zz]
		yy<-yy[zz]
		
	}
	#names(xx)<-names(yy)<-NULL
	if(!ind){xx<-ox[xx];yy<-oref[yy]}
	names(xx)[is.na(xx)]<-nms[is.na(xx)]
	names(yy)[is.na(yy)]<-nms[is.na(yy)]

	if(!identical(names(xx),names(yy)) || !identical(names(xx),nms)){stop("k.indSort: error")}

	list(x=xx,ref=yy,names=nms,ind=ind,inter=inter)
}
indSortAsY<-function(x,y,inter=F,ind=T){

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
k.bindSort<-function(x=NULL,y=NULL,inter=F,opt=c("row","col")){opt<-opt[1]
	if(is_nd1class(x)){x<-as_rowmatrix(x)}
	if(is_nd1class(y)){y<-as_rowmatrix(y)}
	ox<-x;oy<-y
	x<-MakeNames(x);y<-MakeNames(y)
	if(is(x,"data.frame")&&!is(y,"data.frame")){y<-as.data.frame(y)}
	if(!is(x,"data.frame")&&is(y,"data.frame")){x<-as.data.frame(x)}	
	
	if(opt%in%c("row","r")){x<-t(x);y<-t(y)}
	z<-k.indSort(x=colnames(y),ref=colnames(x),inter=inter,ind=T)
	
	nx<-k.add.NA.col.dataframe(x=x,ind=z$ref)#x[,z$ref,drop=F]
	ny<-k.add.NA.col.dataframe(x=y,ind=z$x)#y[,z$x,drop=F]
	stopifnot(ncol(nx)==ncol(ny))
	colnames(nx)<-colnames(ny)<-z$names
	if(opt%in%c("row","r")){nx<-t(nx);ny<-t(ny)}
	list(x=nx,y=ny)
}

k.cbind<-function(x,y=NULL,inter=FALSE,opt=c("row","col")){opt<-opt[1]
	if(length(x)==0){return(y)}
	if(length(y)==0){return(x)}
		
	if(are(list(x,y),c("nd1class","nd2class"))){
		z<-k.bindSort(x=x,y=y,inter=inter,opt=opt)
		if(opt%in%c("row")){zo<-cbind(z$x,z$y)}
		else if(opt%in%c("col")){zo<-rbind(z$x,z$y)}
		else{stop("bad opt")}
		}
	else if(are(list(x,y),'list')){
		x<-MakeNames(x=x,force=FALSE);y<-MakeNames(x=y,force=FALSE)
		z<-k.indSort(x=names(y),ref=names(x),ind=T,inter=T)
		nx<-x[z$ref];ny<-y[z$x]
		stopifnot(identical(names(nx),names(ny)))
		zo<-mapply(x=nx,FUN=k.cbind,y=ny,inter=inter,opt=opt,SIMPLIFY=FALSE)
		}
	return(zo)
}
ncbind<-function(x,y=NULL,inter=FALSE){k.cbind(x=x,y=y,inter=inter,opt="row")}
nrbind<-function(x,y=NULL,inter=FALSE){k.cbind(x=x,y=y,inter=inter,opt="col")}
#
sameXY_names<-function(x=NULL,y=NULL){
	assert.are(list(x,y),c("numeric","logical","character"))
	
	x<-MakeNames(x,nmvar="X")
	y<-MakeNames(y,nmvar="X")
	if(length(y)==0 &&length(x)>0){return(x)}
	if(length(x)==0 &&length(y)>0){return(y)}
	
	if(are(list(x,y),"numeric")){zz<-as.numeric(NA)}
	else if(are(list(x,y),"character")){zz<-as.character(NA)}
	else if(are(list(x,y),"logical")){zz<-NA}

	fint<-as.character(intersect(names(x),names(y)))
	yfint<-as.character(setdiff(names(y),names(x)))
	ml.xy<-length(x)+length(yfint)
	nml.xy<-c(names(x),yfint)

	nx<-rep(zz,ml.xy);names(nx)<-nml.xy;ny<-nx

	nx[names(x)]<-x
	ny[names(y)]<-y
	
	stopifnot(identical(names(nx),names(ny)))
	list(nx=nx,ny=ny)
}
sameAsX_names<-function(x=NULL,y=NULL){
	
	assert.are(list(x,y),c("numeric","logical","character"))
	stopifnot(length(x)>=length(y) && all(names(y)%in%names(x)))
	z<-sameXY_names(x=x,y=y)
	
	z$ny
}
sameAsY<-function(x=NULL,y=NULL){sameAsX_names(x=y,y=x)}
#-------------
nsize<-function(x){
	if(is_any(x,c("list","NULL","numeric","logical","character"))){c(1,length(x))}
	else if(is_any(x,c("matrix","array","data.frame"))){c(nrow(x),ncol(x))}
	else if(is_nothing(x)){c(0,0)}
	else {c(length(x),1)}
}


string2vect<-function (x, sep = ",", fixed  = TRUE) {	
	assert.is(x, "nd1class")
	if(!is(x,'character')){x<-as(x,'chararcter')}
	zo<-strsplit(x=x, split=sep, fixed = fixed , perl = FALSE, useBytes = FALSE)
	sapply(1:length(zo),FUN=function(i){zo[[i]]})
}

 
#-------------- 
nunique<-uniquex<-function(x,y=NULL,vip=1){
	k.nunique<-function(x,y=NULL,vip=1){
		assert.is(vip,'numeric')
		nm<-names(x);unm<-unique(nm)
		if(is_unique(nm)){return(x)}
		#nm.u<-make.unique(nm)
		lnm<-table(nm)
		z<-lapply(1:length(unm),FUN=function(i){
			indx<-which(nm%in%unm[i]);ix<-min(c(vip,length(indx)))
			x[[indx[ix]]]})
		names(z)<-unm;z
        }
	x<-c(x,y)
	if(is_nd1class(x)){zo<-unique(x)}
	else if(is_vide(x)){zo<-x}	
	else if(is(x,'list')){zo<-k.nunique(x=x)}
	#assert.are(list(x),c('list','data.frame','NULL'))
	zo
}

est2list<-function(x){
	
	if (!isS4(x)&&is(x,"list")){zx<-lapply(x,est2list);return(zx)}
	if(!isS4(x)){return(x)}
	slnm<-slotNames(x)
 
	z<-lapply(1:length(slnm),function(i){
	     z1<-eval(parse(text=paste("x@",slnm[i],sep="")))	
	     if(!isS4(z1)){return(z1)}
	     else {z1<-est2list(x=z1)}
	     z1})
	names(z)<-slnm
	assert.is(z,"list")
	if(all(names(z)!=".Data")){
	class(z)<-c("list",paste(class(x),"list",sep="."))}
	z
}

list2est<-function(x,n.object=NULL){
    
	if(is(x,'list') & length(class(x))>1) {
		ind<-grep(x=class(x),pattern='.list',fixed=T)
		if(length(ind)==0){break()}
		obj<-class(x)[ind]
		n.object<-string2vect(obj,sep='.list')[1]}
		
	if(!is(x,'list') || is_vide(n.object)){return(x)}
	nms<-slotNames(n.object)
	assert.is(x,'list')
	if(!all(nms%in%names(x))){message('could not convert list to ',n.object);return(x)}
	nx<-new(n.object);cz<-NULL
	for (i in 1:length(nms)){
		eval(parse(text=paste('cz<-class(nx@',nms[i],')',sep='')))#XXX:no visible binding for global variable "cz"
		if(!isS4(cz)){eval(parse(text=paste('nzi<-list2est(x$',nms[i],',n.object=cz)',sep='')))}
		else{eval(parse(text=paste('nzi<-x$',nms[i],sep='')))}
		eval(parse(text=paste('nx@',nms[i],'<-nzi',sep='')))
  }
  validObject(nx)
  nx
}
#
make_labels<-function(n,nmvar=c("X"),n.ini=1){n.ini<-n.ini[1]
    if(is_vide(n)){return(NULL)}
    stopifnot(n>0)
    paste(nmvar[1],c(1:n)+n.ini-1,sep="")
}

MakeNames<-function(x,nmvar=c("X","I"),force=FALSE,n0=1){
	k.MakeNames<-function(x,col=T,nmvar="I",force=FALSE,unique=T,n0=1){
		if(is_nothing(x)){return(x)}
		k.nms<-function(nm,nn){
			if(is_nothing(nm) | all(is_unk(nm)) | force==T ){
				nm<-make_labels(n=nn,nmvar=nmvar[1],n.ini=n0[1])}
			if(any(is.na(nm)) | any(nchar(nm)==0)){
				ix<-which(is.na(nm) | nchar(nm)==0)
				nm[ix]<-make_labels(n=length(ix),nmvar=nmvar[1],n.ini=n0[1])}
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
		x<-k.MakeNames(x=x,nmvar=nmvar[1],force=force,n0=n0)}
	else if(is_any(x,c("array","matrix","data.frame"))){
		
		x<-k.MakeNames(x,col=T,nmvar=nmvar[max(c(2,length(nmvar)))],force=force,n0=n0)
		x<-k.MakeNames(x,col=F,nmvar=nmvar[1],force=force,n0=n0)
	}
	x
}
list2matrix<-function(x){
	if(is(x,"matrix")){return(x)}
	assert.is(x,"list")
	x<-MakeNames(x)
	nx<-x1<-x[[1]]
	for(i in 2:length(x)){
		nx<-nrbind(nx,x[[i]])
	}
	nx
	
}
matrix2list<-function(x){
	if(is(x,"list")){return(x)}
	if(is_nd1class(x)){x<-as_rowmatrix(x)}
	assert.is(x,"matrix")
	nx<-lapply(1:nrow(x),FUN=function(i){x[i,]})
	names(nx)<-rownames(nx)
	nx
	
}
