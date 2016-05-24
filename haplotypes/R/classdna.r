### class Dna, functions and methods##march 2015, cnrakt#class DnasetClass(Class="Dna", representation=representation(sequence="matrix",seqlengths="numeric",seqnames="character"))#validity: testing Dna objectvalidity <- function(x){	seq<-x@sequence		seq[seq==0]<-"?"	seq[seq==1]<-"A"	seq[seq==2]<-"C"	seq[seq==3]<-"G"	seq[seq==4]<-"T"	seq[seq==5]<-"-"		unchar<- setdiff(c(seq),c("A","C","G","T","a","c","g","t","-","?"))		if(length(unchar)>0)	{				seq[!is.na(match(c(seq),unchar))]<-"?"		if(length(unchar)==1) warning(paste("invalid character ","'", paste(unchar,collapse=","),"'", " was replaced with '?'",sep=""),"\n(valid characters: ", paste(c("A","C","G","T","a","c","g","t","-","?","0","1","2","3","4","5"),collapse=","),")")		if(length(unchar)>1) warning(paste("invalid characters ","'", paste(unchar,collapse=","),"'", " were replaced with '?'",sep=""),"\n(valid characters: ", paste(c("A","C","G","T","a","c","g","t","-","?","0","1","2","3","4","5"),collapse=","),")")	}	x@sequence<-seq	x	}#Set initialize Class DnasetMethod("initialize", "Dna", function(.Object,sequence=matrix(,0,0),seqlengths=integer(0),seqnames=character(0)) {		.Object@sequence<-sequence	.Object@seqlengths<-seqlengths	if(!is.character(seqnames)) seqnames<-as.character(seqnames)	if(nrow(sequence)>0&length(seqnames)==0) seqnames<-as.character(1:nrow(sequence))	if(nrow(sequence)>0&length(seqlengths)==0) .Object@seqlengths<-rep(ncol(sequence),nrow(sequence))	.Object@seqnames<-seqnames	validity(.Object)})#Show method for Dna objectssetMethod("show","Dna", function(object){		cat("*** S4 Object of Class Dna ***\n\n")	cat("\nNumber of DNA sequences: \n")	cat(nrow(object@sequence),"\n")	cat("\nNames: \n")	if(nrow(object@sequence)>6) cat(c(head(object@seqnames),"..."),"\n") else cat(head(object@seqnames),"\n")	cat("\nLength of Shortest DNA Sequence: \n")	if(length(object@seqlengths)==0) cat(integer(0)) else cat(min(object@seqlengths),"\n")	cat("\n\nLength of Longest DNA Sequence: \n")	if(length(object@seqlengths)==0) cat(integer(0)) else cat(max(object@seqlengths),"\n")	cat("\n\nslots of an object Dna:\n")	cat(slotNames(object),"\n")	})#Coerce Dna objects to matrixsetMethod("as.matrix","Dna", function(x){	x@sequence	})#Coerce Dna objects to data.framesetMethod("as.data.frame","Dna", function(x){	df<-as.data.frame(x@sequence)	colnames(df)<-1:ncol(df)	df	})#Coerce Dna objects to listsetMethod("as.list","Dna", function(x) {	n<-nrow(x@sequence)	l<-vector("list",n)	lseq<-x@seqlengths	for(i in 1:n) l[[i]]<-x@sequence[i,1:lseq[i]]	names(l)<-x@seqnames	l	})#Coerce Dna objects to numeric matrixsetMethod("as.numeric","Dna", function(x){	seq<-toupper(x@sequence)	seq[seq=="?"]<-0	seq[seq=="A"]<-1	seq[seq=="C"]<-2	seq[seq=="G"]<-3	seq[seq=="T"]<-4	seq[seq=="-"]<-5	x<-matrix(as.numeric(unlist(seq)),nrow=nrow(seq),ncol=ncol(seq),dimnames=list(x@seqnames,1:ncol(seq)))	x})	#names method for Dna objectssetMethod("names","Dna", function(x) {	x@seqnames		})#names replace method for Dna objectssetReplaceMethod("names","Dna", function(x,value) {	if(is.numeric(value)) value<-as.character(value)	rownames(x@sequence)<-value	x@seqnames<-value	x})#length method for Dna objectssetMethod("length","Dna", function(x) {	ncol(x@sequence)		})#ncol method for Dna objectssetMethod("ncol","Dna", function(x) {	ncol(x@sequence)		})#nrow method for Dna objectssetMethod("nrow","Dna", function(x) {	nrow(x@sequence)		})#Extract method for Dna objectssetMethod("[","Dna", function(x,i=1:nrow(x),j=1:ncol(x),as.matrix=TRUE) {	seq<-x@sequence[i,j,drop=FALSE]	lseq<-x@seqlengths[i]	seqnames<-x@seqnames[i]	if(as.matrix) return(seq) else new("Dna",sequence=seq,seqlengths=lseq,seqnames=seqnames)	})#Extract replace method for Dna objectssetReplaceMethod("[","Dna", function(x,i,j,value) {	x@sequence[i,j]<-value	x<-validity(x)	x})#internal function: read fasta, the first line of the files must start with a ">" (greater-than) symbol.read.fas<-function(file){	fas<-scan(file, what="char",quote="", sep="\n", strip.white=TRUE, quiet= TRUE)		l<-length(fas)	if(!l) stop(paste(file,"file is empty."))	tes<- diff(grep(">",fas))	test1<-any(tes==1)	test2<-any(tes>2)	if(test1&test2) stop("invalid fasta format, description lines of each sequence must start with a '>' (greater-than) symbol and file must not contain empty sequences")	if(test1) stop("invalid fasta format, file must not contain empty sequences")	if(test2) stop("invalid fasta format, description lines of each sequence must start with a '>' (greater-than) symbol")	seqnames<-gsub(">","",fas[seq(1,l,2)])	seqlist<-lapply(fas[seq(2,l,2)],strsplit,"")	seqlist<-lapply(seqlist,unlist)		lseq<-sapply(seqlist,length)	maxseql<-max(lseq)		seq.mat<-matrix("-",l/2,maxseql,dimnames=list(seqnames,1:maxseql))		for(i in 1:(l/2)) seq.mat[i,1:lseq[i]]<-seqlist[[i]]		x<-new("Dna",sequence=seq.mat,seqlengths=lseq,seqnames=seqnames)		x<-validity(x)		x	}		
				

#Generic subssetGeneric (name= "subs",def=function(x,...)standardGeneric("subs"))#subs method for Dna objectssetMethod(f="subs", signature= "Dna", definition=function (x,fifth=FALSE){		seq<-toupper(x@sequence)	if(!fifth) seq[is.na(match(seq,c("A","T","C","G")))]<-"?"	if(fifth) seq[is.na(match(seq,c("A","T","C","G","-")))]<-"?"	uniqs<-apply(seq,2,unique)	whichpoly<-which(unlist(lapply(uniqs,length))>1)	polymat<-seq[,whichpoly,drop=FALSE]		if(ncol(polymat)>0)	{				polylist<- vector("list",0)				for(i in 1:ncol(polymat))		{						substit<-unique(polymat[,i,drop=FALSE][polymat[,i,drop=FALSE]!="?"])						if(length(substit)>1)			{								poly<-list(substit)				names(poly)<-whichpoly[i]				polylist<-c(polylist,poly)			}					}		seq<-seq[,as.numeric(names(polylist)),drop=FALSE]		return(list(subsmat=seq,subs=polylist,subsmnum=length(polylist)))			} else {		return(list(subsmat=polymat,subs=list(),subsmnum=0))	}})# internal function: fillendgapsfillendgaps<-function(x,find="-",replace="?"){	nc<-ncol(x)	for(i in 1:nrow(x))	{		f<-which(x[i,]==find)		fl<-length(f)		if(fl>0) 		{			if(fl==1)			{				beg<-f				end<-f				indell<-1			} else 			{				fpoz<-c(2,f[2:fl]-f[1:(fl-1)])				beg<-f [which(fpoz!=1)]				end<-c(f[which(fpoz[-1]!=1)],f[fl])			}						fl<-length(beg)						if(end[fl]==nc) x[i,beg[fl]:end[fl]]<- replace					}	}		x	}#internal function: alltest

alltest<-function(x,char="-") all(x==char)

# internal function: simple indel coding

sic<-function(x)
{
	
	allindel<-apply(x,2,alltest,char=5)
	x[,allindel]<-0
	
	indelmatrix<-matrix(NA,0,3)
	
	for(i in 1:nrow(x))
	{
		indels<-which(x[i,]==5)
		if(length(indels)>0) 
		{
			if(length(indels)==1)
			{
				beg<-indels
				end<-indels
				indell<-1
			} else 
			{
				indelpoz<-c(2,indels[2:length(indels)]-indels[1:(length(indels)-1)])
				beg<-indels [which(indelpoz!=1)]
				end<-c(indels[which(indelpoz[-1]!=1)],indels[length(indels)])
				indell<-1+end-beg
			}
			
			indelmatrix<-rbind(indelmatrix,cbind(beg,end,indell))
			indelmatrix<-unique(indelmatrix)
			
		}
	}
	
	if(nrow(indelmatrix)==0) return(list(indels=indelmatrix,codematrix=matrix(NA,nrow(x),0)))

	indelmatrix<-indelmatrix[order(indelmatrix[,1],indelmatrix[,2]),,drop=FALSE]
	rownames(indelmatrix)<-1:nrow(indelmatrix)
	
	codematrix<-matrix(0,nrow(x),nrow(indelmatrix))
	element<-vector("list",nrow(indelmatrix))
	
	for(i in 1:ncol(codematrix))
	{
		beg<-indelmatrix[i,1]
		end<-indelmatrix[i,2]
		begother<-indelmatrix[-i,1,drop=TRUE]
		endother<-indelmatrix[-i,2,drop=TRUE]
		element[[i]]<-as.numeric(names(which(beg>=begother&end<=endother)))
		
	}	
	
	
	for(i in 1:ncol(codematrix))
	{
		beg<-indelmatrix[i,1]
		end<-indelmatrix[i,2]
		indelpart<-x[,beg:end]
		if(indelmatrix[i,3]==1) indelpart<-matrix(indelpart,,1)
		indels<- as.numeric(apply(indelpart,1,alltest,char=5))
		
		codematrix[,i]<-indels
		
		el<-element[[i]]
		if(length(el)>0)
		{
			for(j in el)
			{
				beg<-indelmatrix[j,1]
				end<-indelmatrix[j,2]
				indelpart<-x[,beg:end]
				if(indelmatrix[j,3]==1) indelpart<-matrix(indelpart,,1)
				longerindels<- apply(indelpart,1,alltest,char=5)
				codematrix[longerindels,i]<--1
			}
			
		}
		
	}
	
	return(list(indels=indelmatrix,codematrix=codematrix))
	
}
#Generic indelcoder (only simple indel coding is available)setGeneric (name= "indelcoder",def=function(x,...)standardGeneric("indelcoder"))#indelcoder method for Dna objectssetMethod(f="indelcoder", signature= "Dna", definition=function(x){		x<-fillendgaps(x)	seq<-as.numeric(x)	sc<-sic(seq)		return(list(indels=sc$indels,codematrix=sc$codematrix))})#Generic distancesetGeneric (name= "distance",def=function(x,...)standardGeneric("distance"))#distance method for Dna objectssetMethod(f="distance", signature= "Dna", definition=function(x,subset=NULL,indels="sic"){		if(any(subset<1)|any(subset>nrow(x))) stop(paste("elements of 'subset' must be integers in the range [1,",nrow(x),"]",sep=""))		indmet<- c("sic","5th","missing")	matchmet <- pmatch(indels, indmet)	if (is.na(matchmet)) 	stop("invalid indel coding method", paste("", indels))	indels<-indmet[matchmet]		d<-nrow(x)	if(nrow(x)==1) return(dist(0))	dmat<-matrix(NA,d,d)	rownames(dmat)<-x@seqnames			comb<-combn(d,2)		if(!is.null(subset))	{		newcomb<-comb		newcomb[!is.na(match(newcomb,subset))]<-0		comb<-comb[,newcomb[1,]*newcomb[2,]==0]	}		if(indels=="sic")	{				indeldat<-indelcoder(x)$codematrix				if(ncol(indeldat)>0)		{			subsdat<-subs(x,fifth=FALSE)$subsmat						for(i in 1:ncol(comb))			{				ind<-comb[,i]				indeller<-indeldat[ind,,drop=FALSE]				indeller<-indeller[,indeller[1,]!=-1&indeller[2,]!=-1,drop=FALSE]								subslar<-subsdat[ind,,drop=FALSE]				subslar<-subslar[,subslar[1,]!="?"&subslar[2,]!="?",drop=FALSE]								dmat[ind[2],ind[1]]<-sum(indeller[1,]!=indeller[2,])+sum(subslar[1,]!=subslar[2,])							}		} else 		{			indels<-"missing"					}	}			if(indels=="5th")	{		x<-fillendgaps(x)		subsdat<-subs(x,fifth=TRUE)$subsmat		for(i in 1:ncol(comb))		{			ind<-comb[,i]						subslar<-subsdat[ind,,drop=FALSE]			subslar<-subslar[,subslar[1,]!="?"&subslar[2,]!="?",drop=FALSE]						dmat[ind[2],ind[1]]<-sum(subslar[1,]!=subslar[2,])					}					}		if(indels=="missing")	{		subsdat<-subs(x)$subsmat		for(i in 1:ncol(comb))		{			ind<-comb[,i]						subslar<-subsdat[ind,,drop=FALSE]			subslar<-subslar[,subslar[1,]!="?"&subslar[2,]!="?",drop=FALSE]						dmat[ind[2],ind[1]]<-sum(subslar[1,]!=subslar[2,])					}			}			return(as.dist(dmat))})#Generic polymorpsetGeneric (name= "polymorp",def=function(x,...)standardGeneric("polymorp"))#polymorp method for Dna objectssetMethod(f="polymorp", signature= "Dna", definition=function(x,pair,indels="sic"){	if(nrow(x)<2) stop("at least two DNA sequences are required") 	if(length(pair)!=2) stop("'pair' must be vector of length 2")	if(any(pair<1)|any(pair>nrow(x))) stop(paste("elements of 'pair' must be integers in the range [1,",nrow(x),"]",sep=""))		indmet<- c("sic","5th","missing")	matchmet <- pmatch(indels, indmet)	if (is.na(matchmet)) 	stop("invalid indel coding method", paste("", indels))	indels<-indmet[matchmet]		seqmat<-x@sequence	d<-nrow(seqmat)		polylist<-list(list(),list())		names(polylist)<-c("indels","subst")		if(indels=="sic")	{				indel.obj<-indelcoder(x)		indeldat<-indel.obj$codematrix		subs.obj<-subs(x)		subsdat<-subs.obj$subsmat				indeller<-indeldat[pair,,drop=FALSE]				nonmis<-which(indeller[1,]!=-1&indeller[2,]!=-1)		indeller<-indeller[,nonmis,drop=FALSE]		diff<-which(indeller[1,]!=indeller[2,])		polyind<-nonmis[diff]		ins<-indel.obj$indels[polyind,,drop=FALSE]		indellist<-vector("list",0)		if(nrow(ins)>0)		{			indellist<-vector("list",nrow(ins))			names(indellist)<-ins[,1]			for(i in 1:nrow(ins)) indellist[[i]]<-  seqmat[pair,ins[i,1]:ins[i,2],drop=FALSE]			polylist[[1]]<-indellist		}				subslar<-subsdat[pair,,drop=FALSE]		polysubs<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!=subslar[2,],drop=FALSE]			subsmat<-polysubs		subslist<-vector("list",0)		if(ncol(subsmat)>0)		{			subslist<-vector("list",ncol(subsmat))			names(subslist)<-colnames(subsmat)			for(i in 1:ncol(subsmat)) subslist[[i]]<- subsmat[,i,drop=FALSE]			}		polylist[[2]]<-subslist					}			if(indels=="5th")	{				x<-fillendgaps(x)		subs.obj<-subs(x,fifth=TRUE)		subsdat<-subs.obj$subsmat		subslar<-subsdat[pair,,drop=FALSE]		polysubs<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!=subslar[2,],drop=FALSE]		subsmat<-polysubs				subslist<-vector("list",0)		if(ncol(subsmat)>0)		{						subslist<-vector("list",ncol(subsmat))			names(subslist)<-colnames(subsmat)			for(i in 1:ncol(subsmat)) subslist[[i]]<- subsmat[,i,drop=FALSE]			}			polylist[[2]]<-subslist	}		if(indels=="missing")	{		subs.obj<-subs(x)		subsdat<-subs.obj$subsmat		subslar<-subsdat[pair,,drop=FALSE]		subslar<-subsdat[pair,,drop=FALSE]		polysubs<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!=subslar[2,],drop=FALSE]				subsmat<-polysubs			subslist<-vector("list",0)		if(ncol(subsmat)>0)		{						subslist<-vector("list",ncol(subsmat))			names(subslist)<-colnames(subsmat)			for(i in 1:ncol(subsmat)) subslist[[i]]<- subsmat[,i,drop=FALSE]			}			polylist[[2]]<-subslist	}				return(polylist)})#Generic as.dnasetGeneric (name= "as.dna",def=function(x,...)standardGeneric("as.dna"))#Coerce Haplotype objects to Dna objectsetMethod(f="as.dna", signature= "Haplotype", definition=function(x){		seq<-x@sequence	if(!nrow(seq)) stop("Haplotype object does not contain DNA sequences")	dnaobj<-new("Dna",sequence=seq,seqlengths=ncol(seq),seqnames=rownames(seq))	dnaobj})#Coerce matrix to Dna objectsetMethod(f="as.dna", signature= "matrix", definition=function(x){	if(is.numeric(x))	{		x[x==0]<-"?"		x[x==1]<-"A"		x[x==2]<-"C"		x[x==3]<-"G"		x[x==4]<-"T"		x[x==5]<-"-"		x[x==NA]<-NA	}		dnaobj<-new("Dna",sequence=x,seqlengths=rep(ncol(x),nrow(x)),seqnames=rownames(x))	dnaobj})#Coerce data.frame to Dna objectsetMethod(f="as.dna", signature= "data.frame", definition=function(x){	if(is.numeric(x))	{		x[x==0]<-"?"		x[x==1]<-"A"		x[x==2]<-"C"		x[x==3]<-"G"		x[x==4]<-"T"		x[x==5]<-"-"		x[x==NA]<-NA	}		dnaobj<-new("Dna",sequence=as.matrix(x),seqlengths=rep(ncol(x),nrow(x)),seqnames=rownames(x))	dnaobj})#Coerce list to Dna objectsetMethod(f="as.dna", signature= "list", definition=function(x){	lseq<-sapply(x,length)	l<-length(x)	seqnames<-names(x)	maxseql<-max(lseq)		seq.mat<-matrix("-",l,maxseql,dimnames=list(seqnames,1:maxseql))	for(i in 1:(l)) seq.mat[i,1:lseq[i]]<-x[[i]]		dnaobj<-new("Dna",sequence=seq.mat,seqlengths=lseq,seqnames=seqnames)	dnaobj})#append method for Dna objects setMethod(f="append", signature= "Dna", definition=function(x,values){	seq<- rbind(x@sequence,values@sequence)	dnaobj<-new("Dna",sequence=seq,seqlengths=rep(ncol(x),nrow(x)),seqnames=c(x@seqnames,values@seqnames))	dnaobj})# internal function: pairwise nei raw D (D) : Nei s average number of differences between populations  (Nei and Li, 1979)pairneidist<-function(x,populations){		x<-as.matrix(x)			pops<-unique(populations)	npop<-length(pops) 	withinnei<-vector("numeric",npop)			for(i in 1:npop)	{		p<-populations==pops[i]		xw<-x[p,p]		n<-ncol(xw)		if(!is.null(ncol(xw)))  nei <-mean(as.dist(xw)) else nei <- 0 		withinnei[i]<-nei	}		neidistmat<-matrix(0,npop,npop)		comb<-combn(npop,2)		for(i in 1:ncol(comb))	{				p1<-populations==pops[comb[1,i]]		p2<-populations==pops[comb[2,i]]		nd<-mean(x=as.matrix(x[p1,p2]))		neidistmat[comb[2,i],comb[1,i]]<-nd	}		listele<-list(Between=neidistmat,Within=withinnei,Uniquepopulations=pops)		return(listele)	}#Generic pairneisetGeneric (name= "pairnei",def=function(x,...)standardGeneric("pairnei"))#pairnei method for Dna objectsetMethod("pairnei","Dna", function(x,populations,indels="sic"){		d<-distance(x,indels=indels)	if(length(d)==0) stop("at least two DNA sequences are required")	if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")	if(length(populations)>nrow(x)) stop(paste("'populations' must be vector of length",nrow(x),"(equal or less than the number of DNA sequences)",sep=" "))	n<-pairneidist(d,populations)	n})#pairnei method for matrix objectsetMethod("pairnei","matrix", function(x,populations){		if(nrow(x)==1) stop("dimensions of the distance matrix (x) must be greater than one")	if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")	if(length(populations)>nrow(x)) stop(paste("'populations' must be vector of length",nrow(x),"(equal or less than the dimensions of the distance matrix 'x')",sep=" "))	n<-pairneidist(as.dist(x),populations)	n})#pairnei method for dist objectsetMethod("pairnei","dist", function(x,populations){	if(length(x)==0) stop("length of the distance object (x) must be greater than zero")	if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")	if(length(populations)>nrow(as.matrix(x))) stop(paste("'populations' must be vector of length",nrow(as.matrix(x)),"(equal or less than the dimensions of the distance matrix 'as.matrix(x)')",sep=" "))	n<-pairneidist(x,populations)	n	})