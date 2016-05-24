#Utility to do set analysis
#Author: Minghui Wang, minghui.wang@mssm.edu
#Date: 20 July, 2014

#
jaccard=function(x){
	if(!is.list(x)) stop('Input x must be list\n')
	nL=length(x)
	if(nL<2) stop('Input x should have at least two entries\n')
	x=lapply(x,unique)
	Mat=matrix(NA,nL,nL)
	colnames(Mat)=rownames(Mat)=names(x)
	diag(Mat)=1
	for(i in 1:(nL-1)){
		for(j in (i+1):nL) Mat[i,j]=Mat[j,i]=sum(x[[i]] %in% x[[j]])/length(c(x[[i]],x[[j]]))
	}
	Mat
}
#
enumerateIntersecSizes=function(x,degree=NULL){
	if(!is.null(degree)) return(incIntersect(x,degree))
	otab=exclusiveIntersect0(x)
	exc2incIntersect(otab)
}
#list all possible intersections
intersectElements=function(x){
#return Venn diagram entry sizes
#x: a list of sets
	if(!is.list(x)) stop('Input x must be list\n')
	nL=length(x)
	if(nL<2) stop('Input x should have at least two entries\n')
	allE=unique(unlist(x))
	barcodes=rep('',length(allE))
	for(i in 1:nL){
		barcodes=paste(barcodes,ifelse(allE %in% x[[i]],'1','0'),sep='')
	}
	data.frame(Entry=allE,barcode=barcodes,stringsAsFactors=FALSE)
}
#enumerate all overlap and non-overlap set sizes
exclusiveIntersect0=function(x){
#x: a list of sets
#return mutual exclusive intersection sizes including empty intersections
	intersects=intersectElements(x)
	nL=length(x)
	barcodes=mkBarcode(nL)
	otab=sapply(barcodes,function(a) 0)
	tab=table(intersects$barcode)
	otab[names(tab)]=tab
	otab
}
exc2incIntersect=function(x){
#x, an object generated from function exclusiveIntersect0
#return inclusive subset sizes
	otab=x
	otab[]=0
	C1 = lapply(strsplit(names(x),''), function(c11) c11 == '1')
	for(i in 1:length(x)){
		a=C1[[i]]
		rel=sapply(C1,function(b) all(a[b]==TRUE))
		otab[rel]=otab[rel]+x[i]
	}
	otab
}
#reverse barcode
deBarcode <- function(barcode,grp){
	sapply(barcode,function(b){
		s=grp[strsplit(b,'')[[1]] == '1']
		s=paste(s,collapse=' & ')
		s
	})
}
#compute intersection sizes for given overlap degree
incIntersect=function(x,degree=NULL){
#x is a list of sets
	if(!is.list(x)) stop('Input x must be list\n')
	nL=length(x)
	if(nL<2) stop('Input x should have at least two entries\n')
	allE=unique(unlist(x))
	nE=length(allE)
	BinMat=sapply(x, function(d) allE %in% d )
	barcodes=mkBarcode.degree(nL,degree)
	otab=sapply(barcodes,function(a) 0)
	for(i in 1:length(otab)){
		i1=which(strsplit(names(otab)[i],'')[[1]] == '1')
		otab[i]=sum(rowSums(BinMat[,i1,drop=FALSE]) == length(i1))
	}
	otab
}
#
intersect=function(x,y,...){
	dat=list(x,y,...)
	if(length(dat)<2) return(unlist(dat))
	common=as.vector(dat[[1]])
    for(i in 2:length(dat)){
		common=unique(common[match(as.vector(dat[[i]]), common, 0L)])
		if(length(common)==0) break
	}
	common
}
union=function(x,y,...){
	dat=list(x,y,...)
	if(length(dat)<2) return(unlist(dat))
	u=as.vector(dat[[1]])
    for(i in 2:length(dat)){
		u=unique(c(u,as.vector(dat[[i]])))
	}
	u
}
intersect.list=function(x){
	if(! is.list(x)) stop('Input must be a list\n')
	if(length(x)<2) return(unlist(x))
	common=as.vector(x[[1]])
    for(i in 2:length(x)){
		common=unique(common[match(as.vector(x[[i]]), common, 0L)])
		if(length(common)==0) break
	}
	common
}
union.list=function(x){
	if(! is.list(x)) stop('Input must be a list\n')
	if(length(x)<2) return(unlist(x))
	u=as.vector(x[[1]])
    for(i in 2:length(x)){
		u=unique(c(u,as.vector(x[[i]])))
	}
	u
}
