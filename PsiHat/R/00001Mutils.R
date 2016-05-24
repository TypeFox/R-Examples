
setGeneric("CorrectLimits", function(x,...) standardGeneric("CorrectLimits"))
setMethod("CorrectLimits", signature(x = "numeric" ), function(x,maxlim=1,minlim=0){
	 stopifnot(!are_null(list(maxlim,minlim)))
	x[x<minlim]<-minlim
	x[x>maxlim]<-maxlim
	x})
#
vect2string<-function(x,sep="",...){paste(x,sep="",collapse=sep,...)}
string2char<-function(x,...){strsplit(x,split=NULL,...)}
grep_or<-function(x,pattern,fixed=FALSE,exact=FALSE,ind=T,unik=T,...){
	z<-lapply(1:length(pattern),FUN=function(i){
		if(exact){indx<-which(x%in%pattern[i])}
		else{indx<-grep(x=x,pattern=pattern[i],fixed=fixed,...)}
		if(!ind){x[indx]} else {indx}})
	zo<-unlist(z)
	if(unik){zo<-unique(zo)}
		zo
}
grepl_or<-function(x,pattern,fixed=FALSE,exact=FALSE,unik=T,...){
	z<-grep_or(x=x,pattern=pattern,fixed=fixed,exact=exact,ind=T,unik=T,...);length(z)>0}
#
###---------
undefAsNA<-function(x){
	if(is(x,'numeric')){ind<-which(!is.finite(x));if(length(ind)>0){x[ind]<-NA} else {x}}
	else if(is(x,'matrix')){x<-apply(x,2,undefAsNA)}
	x
}
removeRC.from.matrix<-function(x,opt='NA',indx=FALSE){
	ox<-x;if(is_vide(x)){message('empty matrix');return(NULL)}
	if (is_nd1class(x)){
		x<-as_rowmatrix(x)}
		all.eq<-function(x){if(length(table(x))==1){1} else{0}
	}
	
	all_na<-function(x){if(length(which(is.na(x)))==length(x)){1} else{0}}
	all_nonfinit<-function(x){if(length(which(!is.finite(x)))==length(x)){1} else{0}}
	zr.fun<-function(x,fun){
		assert.is(x,"nd2class")
		zr<-vapply(1:nrow(x),FUN.VALUE=numeric(1),FUN=function(i){fun(x[i,])})
		names(zr)<-rownames(zr);zr
	}
	zc.fun<-function(x,fun){assert.is(x,"nd2class")
		zr<-vapply(1:ncol(x),FUN.VALUE=numeric(1),FUN=function(i){fun(x[,i])})
		names(zr)<-colnames(zr);zr
	}
	rem_row<-function(x,zr){if(1 %in% zr){x<-x[-which(zr==1),,drop=FALSE]}
		x}
	rem_col<-function(x,zc){if(1 %in% zc){x<-x[,-which(zc==1),drop=FALSE]}
		x}

	if(opt%in%c('NA','na'))	{zr<-zr.fun(x=x,fun=all_na);zc<-zc.fun(x=x,fun=all_na)}
	else if(opt%in%c('eq','eqs','EQ','EQS')){zr<-zr.fun(x=x,fun=all.eq);zc<-zc.fun(x=x,fun=all.eq)}
	else if(opt%in%c('nfinit','Inf','inf','NaN','nan')){zr<-zr.fun(x=x,fun=all_nonfinit);zc<-zc.fun(x=x,fun=all_nonfinit)}
	else{stop('bad opt')}
	
	x<-rem_row(x=x,zr=zr);indr<-which(zr==1)
	x<-rem_col(x=x,zc=zc);indc<-which(zc==1)
	if(is_nd1class(ox)){nms<-colnames(x);x<-as.numeric(x);names(x)<-nms}
	message(paste('removing',length(which(zr==1)),'rows and',length(which(zc==1)),'columns\n') )
	if(indx){list(matrix=x,indNA.rows=indr,indNA.cols=indc)}
	else{x}
}
#
removeNA.from.matrix<-function(x,indx=FALSE){removeRC.from.matrix(x=x,opt='nfinit',indx=indx)}#x<-undefAsNA(x);
removeEQ.from.matrix<-function(x,indx=FALSE){removeRC.from.matrix(x=x,opt='eq',indx=indx)}
#removeNF.from.matrix<-function(x,indx=FALSE){removeRC.from.matrix(x=x,opt='nfinit',indx=indx)}
#
prep.2matrices<-function(x,y=NULL,paired=FALSE,rm.na=T){
	if(length(y)==0){zo<-list(x=x,y=NULL,info=list(paired=FALSE));return(zo)}
	if(is_nd1class(x)){x<-as_rowmatrix(x)}
	if(is_nd1class(y)){y<-as_rowmatrix(y)}
	assert.are(list(x,y),'matrix')
	nmvarx<-c('X','I')
	if(paired){nmvary<-c('X','I')} else {nmvary<-c('X','S')}
	x<-MakeNames(x,nmvar=nmvarx)                 
	if(rm.na){x<-removeNA.from.matrix(x=x,indx=FALSE)}
		
	if(!is_vide(y)){
		y<-MakeNames(y,nmvar=nmvary)
		if(rm.na){y<-removeNA.from.matrix(x=y,indx=FALSE)}
		indx<-k.indSort(x=rownames(y),ref=rownames(x),ind=T,inter=T)#,rm.na=T,inter=T
		if(length(indx)==0){stop('different rownames (feature names) in x and y')}
		y<-y[indx$x,,drop=FALSE];x<-x[indx$ref,,drop=FALSE]
		if(setequal(colnames(y),colnames(x))||paired==TRUE){
			indx<-k.indSort(x=colnames(y),ref=colnames(x),ind=T,inter=T)
			if(is_vide(indx)){stop('no paired x and y')}
			y<-y[,indx$x,drop=FALSE];x<-x[,indx$ref,drop=FALSE]
			stopifnot(identical(colnames(x),colnames(y)))
		}
		}
	stopifnot(identical(rownames(x),rownames(y)))
	list(x=x,y=y,info=list(paired=paired))
}
setGeneric("is_paired", function(x,y) standardGeneric("is_paired"))
setMethod("is_paired", signature(x = "nd2class",y = "nd2class"), function(x,y){
    x<-MakeNames(x)
    y<-MakeNames(y)
    if(identical(rownames(x),rownames(y)) && identical(colnames(x),colnames(y))){return(TRUE)}
    else{return(FALSE)}
    })
setMethod("is_paired", signature(x = "nd1class",y = "nd1class"), function(x,y){
    x<-as_rowmatrix(x)
    y<-as_rowmatrix(y)
    is_paired(x=x,y=y)})
setMethod("is_paired", signature(x = "nd2class",y = "nd1class"), function(x,y){y<-as_rowmatrix(y);is_paired(x=x,y=y)})
setMethod("is_paired", signature(x = "nd1class",y = "nd2class"), function(x,y){x<-as_rowmatrix(x);is_paired(x=x,y=y)})

#
k.get.from.funtest<-function(x,y=NULL,test.fun=t.test,opt=c("p.value","statistic","parameter","conf.int"),rm.na=TRUE,...){
	fun<-test.fun
	x<-as.numeric(x);if(!is_vide(y)){y<-as.numeric(y)}
	assert.are(list(x,y),c('numeric','NULL'))
	if(any(opt%in%c("p.value","pval","pvalue","pv"))){opt[opt%in%c("p.value","pval","pvalue","pv")]<-"p.value"}
	if(any(opt%in%c("statistic","stat","st","stats"))){opt[opt%in%c("statistic","stat","st","stats")]<-"statistic"}

	if(rm.na){x<-x[!is.na(x)];y<-y[!is.na(x)]}
	if (length(opt)==1){	
		if(is_vide(y)){eval(parse(text=paste('zo<-try(fun(x=x,...)$',opt,')',sep='')))}
		else{eval(parse(text=paste('zo<-try(fun(x=x,y=y,...)$',opt,')',sep='')))}#y<-y[!is.na(y)];
		if(is_err(zo)){zo<-as.numeric(NA);names(zo)<-opt}
	}
	else{
		if(is_vide(y)){eval(parse(text=paste('zo<-try(fun(x=x,...))',sep='')))}
		else{eval(parse(text=paste('zo<-try(fun(x=x,y=y,...))',sep='')))}#y<-y[!is.na(y)];
		if(is_err(zo)){zo<-as.numeric(rep(NA,length(opt)));names(zo)<-opt}
		else{
			if(!is(zo,"list")){zo<-as.list(zo)}
			zo<-zo[opt]
			nzo<-sapply(1:length(zo),FUN=function(i){
				nopt<-rep(opt[i],length(zo[[i]]))
				names(zo[[i]])<-paste(names(zo[[i]]),nopt,sep="");zo[[i]]})
			zo<-nzo
			if(is(zo,"list")){
				zo<-do.call(c,zo)
				names(zo)<-gsub(x=names(zo),pattern="x.",replacement="",fixed =T)
				names(zo)<-gsub(x=names(zo),pattern="y.",replacement="",fixed =T)}
		}
	}
			zo
}

k.get.from.funtest.mat<-function(x,y=NULL,test.fun=t.test,paired=FALSE,opt='p.value',...){rm.na<-FALSE
	xx<-x;yy<-y
	k.i.row<-function(i){
		xx<-x[i,];yy<-y[i,]
		if(paired && identical(names(xx),names(yy))){
			ix<-intersect(which(!is.na(xx)),which(!is.na(yy)))
			if(length(ix)==0){warning('non finite elements');return(as.numeric(NA))}
			xx<-xx[ix];yy<-yy[ix]
			}
		
		do.call(k.get.from.funtest,c(list(x=xx,y=yy,test.fun=test.fun,opt=opt,rm.na=rm.na),argk))
		}
	#-------------
	if(is_any(x,c("nxprnSet"))){xx<-nexprs(x)}
	if(is_any(y,c("nxprnSet"))){yy<-nexprs(y)}
	if(is_any(x,c("nxprnSetPair"))){xx<-nexprs(x);yy<-nexprs(y)}
    	x<-xx;y<-yy
	if(is_nd1class(x)){x<-as_rowmatrix(x)}
	if(is_nd1class(y)){y<-as_rowmatrix(y)}
	
	assert.is(x,'matrix');assert.is(y,c('matrix','NULL'))
	x<-MakeNames(x);y<-MakeNames(y)
	argk<-nunique(list(...))
	
	zz<-prep.2matrices(x=x,y=y,paired=paired,rm.na=T);x<-zz$x;y<-zz$y

	zo<-NULL
	pp1<-k.i.row(1)
	npp1<-length(pp1)
	na.pp1<-rep(as.numeric(NA),npp1);names(na.pp1)<-names(pp1)
	if(nrow(x)>1){
		zo<-vapply(2:nrow(x),FUN.VALUE=numeric(npp1),FUN=function(i){
			zaux<-k.i.row(i)
			if(length(zaux)<npp1){message("problem in row ",i);zaux<-na.pp1}
			zaux})}
	if(length(pp1)==1){zo<-c(pp1,zo)}
	else{zo<-cbind(pp1,zo)}
	if(is(zo,"numeric")){names(zo)<-rownames(x)}
	else if(is(zo,"matrix")){colnames(zo)<-rownames(x)}
	
zo
}
#
get_other.from.testfun<-function(x,y=NULL,test.fun=t.test,paired=FALSE,opt='parameter',...){
	k.get.from.funtest.mat(opt=opt,x=x,y=y,test.fun=test.fun,paired=paired,...)
}
get_pvalues<-function(x,y=NULL,test.fun=t.test,paired=FALSE,...){
	zo<-k.get.from.funtest.mat(opt='p.value',x=x,y=y,test.fun=test.fun,paired=paired,...)
	stopifnot(all(is_prob(zo[is.finite(zo)])))
	zo
}
get_stats<-function(x,y=NULL,test.fun=t.test,paired=FALSE,...){
	k.get.from.funtest.mat(opt='statistic',x=x,y=y,test.fun=test.fun,paired=paired,...)
}

stat2pval<-function(x,pFUN=pt,alternative='two.sided',sym.distrib=T,...){#alternative = "greater" is the alternative that x has a larger mean than y
	x<-MakeNames(x)
	if(sym.distrib&&identical(pFUN,'pchisq')&&alternative%in%c('t','two.sided')){warning('setting sym.distrib=FALSE');sym.distrib<-FALSE}
	if(!sym.distrib&&identical(pFUN,'pt')&&alternative%in%c('t','two.sided')){warning('setting sym.distrib=TRUE');sym.distrib<-T}
	
	if(alternative%in%c('g','greater','t','two.sided')){lower.tail<-FALSE}
	else if(alternative%in%c('l','less')){lower.tail<-TRUE}
	else{stop('bad alternative')}
	k.stat2pval<-function(y,...){
		assert.is(y, "numeric")
		pv<-pFUN(q=y,lower.tail=lower.tail,...)
		if(alternative%in%c('t','two.sided')&&sym.distrib){pv<-2*pmin(pv,pFUN(q=y,lower.tail=!lower.tail,...))}
		else if(alternative%in%c('t','two.sided')&&!sym.distrib){
			pv<-pFUN(q=y,lower.tail=FALSE,...)+pFUN(q=-y,lower.tail=T,...)}
			pv<-pmin(pv,1)
	
		stopifnot(all(is_prob(pv[is.finite(pv)])))
		names(pv)<-names(y)
		pv}
	#--------------------------------
	if(is_nd1class(x)){zo<-k.stat2pval(y=x,...);return(zo)}
	else if(is_nd2class(x)){
		zo<-vapply(1:nrow(x),FUN.VALUE=numeric(ncol(x)),FUN=function(i){
			k.stat2pval(y=x[i,],...)})
		rownames(zo)<-colnames(x)
		colnames(zo)<-rownames(x)
		t(zo)
		
	}
}
pval2stat<-function(x,qFUN=qt,alternative='greater',...){
	x<-MakeNames(x)
		
	if(alternative%in%c('g','greater','t','two.sided')){lower.tail<-FALSE}
	else if(alternative%in%c('l','less')){lower.tail<-TRUE}
	else{stop('bad alternative')}
	if(alternative%in%c('t','two.sided')){warning('wrong statistics from two.sided p-values');x<-x/2}

	k.pval2stat<-function(y,...){
		assert.is(y, "numeric")
		stopifnot(all(is_prob(y[is.finite(y)])))
			
		st<-qFUN(p=y,lower.tail=lower.tail,...)
		names(st)<-names(y)
		st}
	#--------------------------------
	if(is_nd1class(x)){zo<-k.pval2stat(y=x,...);return(zo)}
	else if(is_nd2class(x)){
		zo<-vapply(1:nrow(x),FUN.VALUE=numeric(ncol(x)),FUN=function(i){
			k.pval2stat(y=x[i,],...)})
		rownames(zo)<-colnames(x)
		colnames(zo)<-rownames(x)
		t(zo)	
	}
}
###-------------------------
