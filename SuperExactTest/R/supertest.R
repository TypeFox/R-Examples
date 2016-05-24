#Compute overlap test and visualize intersections between multiple sets
#Author: Minghui Wang
#minghui.wang@mssm.edu
#
setGeneric("supertest", function(x, n=NULL,degree=NULL,...)
           standardGeneric("supertest"),signature="x")
supertest.list<-function(x,n=NULL,degree=NULL,...){
#x, a list of sets
#n, integer, background population size
#..., additional arguments, not implemented

	if(!is.list(x)) stop('Input must be a list\n')
	if(is.null(names(x))) names(x)=paste('Set',1:length(x),sep='')
	set.names=names(x)
	if(any(set.names=='')) stop('Please specify names for each list entry\n')
	obj=list()
	obj$x=x
	obj$set.names=set.names
	obj$set.sizes=sapply(x,function(x) length(unique(x)))
	obj$n=n
	obj$overlap.sizes=enumerateIntersecSizes(x,degree=degree)
	obj$P.value=NULL
	if(!is.null(n)){
		if(any(obj$set.sizes>n)) stop('Background population size should not be smaller than set size\n')
		obj$P.value=sapply(1:length(obj$overlap.sizes),function(i){
			which.set=which(strsplit(names(obj$overlap.sizes)[i],'')[[1]]=='1')
			if(length(which.set)==1) return(NA)
			if(obj$overlap.sizes[i]==0) return(1)
			cpsets(max(obj$overlap.sizes[i]-1,0),obj$set.sizes[which.set],n,lower.tail=FALSE)
		})
		names(obj$P.value)=names(obj$overlap.sizes)
	}
	class(obj)='msets'
	obj
}
setMethod("supertest", signature=c(x="list"), supertest.list)

print.msets=function(x,...){
	cat('A msets object\n')
}
summary.msets=function(object, degree=NULL, ...){
	nL=length(object$x)
	otab=object$overlap.sizes
	if(is.null(degree)) degree=1:nL
	if(any(degree < 1) || any(degree > nL)) stop('Invalid degree value\n')
	odegree=sapply(names(otab),function(d) countCharOccurrences('1',d))
	otab=otab[odegree %in% degree]
	if(length(otab)==0) stop('No data for output\n')
	Barcode=names(otab)
	odegree=odegree[Barcode]
	etab=rep(NA,length(otab))
	if(!is.null(object$n)){
		for(i in 1:length(otab)){
			if(odegree[i] == 1) next
			s=strsplit(Barcode[i],'')[[1]] == '1'
			#if(sum(s)==1) next
			etab[i]=object$n*do.call('prod',as.list(object$set.sizes[s]/object$n))
		}
	}
	res=list(Barcode=Barcode,otab=otab,etab=etab,set.names=object$set.names,set.sizes=object$set.sizes,n=object$n,P.value=object$P.value)
	#find intersections
	el=intersectElements(object$x)
	bc=strsplit(el$barcode,'')
	Elements=sapply(Barcode,function(d){
		id=which(strsplit(d,'')[[1]]=='1')
		od=sapply(bc,function(b) ifelse(all(b[id]=='1'),TRUE,FALSE))
		paste(el[od,1],collapse=', ')
	})
	if(is.null(object$n)){
		res$Table=data.frame(Intersections=deBarcode(Barcode,object$set.names),Degree=odegree,Observed.Overlap=otab,Elements=Elements,stringsAsFactors=FALSE)
	}else{
		res$Table=data.frame(Intersections=deBarcode(Barcode,object$set.names),Degree=odegree,Observed.Overlap=otab,Expected.Overlap=etab,FE=otab/etab,P.value=object$P.value[Barcode],Elements=Elements,stringsAsFactors=FALSE)
	}
	rownames(res$Table)=Barcode
	class(res)='summary.msets'
	res
}
print.summary.msets=function(x,...){
	cat('A msets object with',length(x$set.names),'sets:',x$set.names,'\n')
	if(!is.null(x$n)) cat('Background size:',x$n,'\n')
	cat('Summary of intersections:\n')
	x$Table$Elements=sapply(x$Table$Elements,function(d){
		if(nchar(d)>20) d=paste(substr(d,0,20),' ...',sep='')
		d
	})
	print(x$Table)
}
#
