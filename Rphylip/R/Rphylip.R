## internal function read.multi.dna
## written by Liam J. Revell 2014

read.multi.dna<-function(file,N,n){
	FF<-readLines(file)
	skip<-grep(pattern=paste("   ",n,sep=""),FF)-1
	X<-lapply(skip,read.dna,file=file,format="sequential")
	return(X)
}

## internal function read.multi.phylip.data
## written by Liam J. Revell 2014

read.multi.phylip.data<-function(file,N,n){
	FF<-readLines(file)
	skip<-grep(pattern=paste("   ",n,sep=""),FF)-1
	X<-lapply(skip,read.phylip.data,file=file,format="sequential")
	return(X)
}

## internal function read.multi.rest.data
## written by Liam J. Revell 2014

read.multi.rest.data<-function(file,N,n){
	FF<-readLines(file)
	skip<-grep(pattern=paste("   ",n,sep=""),FF)-1
	X<-lapply(skip,read.rest.data,file=file)
	return(X)
}

## read.phylip.data
## written by Liam J. Revell 2014

read.phylip.data<-function(file,format="interleaved",skip=0,nlines=0,comment.char="#",as.character=FALSE){
	X<-read.dna(file,format,skip,nlines,comment.char,as.character=TRUE)
	class(X)<-"phylip.data"
	X
}

## read.rest.data
## written by Liam J. Revell 2014

read.rest.data<-function(file,skip=0){
	X<-readLines(file)
	nn<-strsplit(X[skip+1]," ")[[1]]
	nn<-as.numeric(nn[nn!=""])
	X<-X[1:(nn[1]*ceiling(nn[2]/50))+1]
	X<-sapply(1:nn[1],function(i,x,nn) paste(x[1:ceiling(nn[2]/50)+(i-1)*ceiling(nn[2]/50)],collapse=""),x=X,nn=nn)
	X<-strsplit(X,"")
	labels<-sapply(lapply(X,function(x) x[1:10]),function(x) x[1:max(setdiff(1:10,which(x==" ")))])
	X<-lapply(lapply(X,function(x) x[11:length(x)]),function(x) x<-x[x!=" "])
	names(X)<-labels
	attr(X,"nenzymes")<-nn[3]
	attr(X,"nsites")<-nn[2]
	class(X)<-"rest.data"
	X
}

## calls seqboot from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rseqboot<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("seqboot")
	if(is.null(path)) stop("No path provided and was not able to find path to seqboot")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet){
		files<-c("ancestors","categories","infile","mixture","outancestors","outcategories","outfile",
			"outmixture","weights")
		if(file.warn(files)==0) return(NULL)
	}
	oo<-vector()
	if(hasArg(type)){ 
		type<-list(...)$type
		type<-tolower(type)
	} else {
		if(class(X)=="DNAbin"||class(X)=="proseq") type<-"sequence"
		else if(class(X)=="phylip.data") type<-"morph"
		else if(class(X)=="rest.data") type="rest"
		else if(class(X)=="matrix") type="gene.freq"
	}
	if(type=="morph") oo<-c(oo,"d")
	else if(type=="rest") oo<-c(oo,rep("d",2),"e")
	else if(type=="gene.freq") oo<-c(oo,rep("d",3))
	if(hasArg(method)) method<-list(...)$method
	else method<-"bootstrap"
	method<-tolower(method)
	if(method=="jacknife") oo<-c(oo,"j")
	else if(method=="permute") oo<-c(oo,rep("j",2))
	if(hasArg(percentage)) percentage<-list(...)$percentage
	else percentage<-100
	if(percentage!=100) oo<-c(oo,"%",percentage)
	if(hasArg(block.size)) block.size<-list(...)$block.size
	else block.size<-1
	if(block.size!=1) oo<-c(oo,"b",block.size)
	if(hasArg(replicates)) replicates<-list(...)$replicates
	else replicates<-100
	if(replicates!=100) oo<-c(oo,"r",replicates)
	if(hasArg(weights)){
		weights<-list(...)$weights
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(hasArg(rate.categories)&&type=="sequence"){
		rate.categories<-list(...)$rate.categories
		write(paste(rate.categories,collapse=""),file="categories")
		oo<-c(oo,"c")
	} else rate.categories<-NULL
	if(hasArg(mixture)&&type=="morph"){
		mixture<-list(...)$mixture
		oo<-c(oo,"x")
		mixture<-toupper(mixture)
		write(paste(mixture,collapse=""),file="mixture")
	} else mixture<-NULL
	if(hasArg(ancestors)&&type=="morph"){
		ancestors<-list(...)$ancestors
		oo<-c(oo,"n")
		ancestors<-toupper(ancestors)
		write(paste(ancestors,collapse=""),file="ancestors")
	} else ancestors<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"i","y",sample(seq(1,99999,by=2),1),"r")
	if(type=="sequence"||type=="morph") write.dna(X)
	else if(type=="rest") write.rest.data(X)
	else if(type=="gene.freq") write.continuous(X)
	system("touch outfile")
	if(!is.null(ancestors)){ 
		system("touch outancestors")
		oo<-c(oo,"r")
	}
	if(!is.null(mixture)){ 
		system("touch outmixture")
		oo<-c(oo,"r")
	}
	if(!is.null(rate.categories)){
		system("touch outcategories")
		oo<-c(oo,"r")
	}
	system(paste(path,"/seqboot",sep=""),input=oo,show.output.on.console=(!quiet))
	if(type=="sequence") XX<-read.multi.dna(file="outfile",N=replicates,n=nrow(X))
	if(type=="morph") XX<-read.multi.phylip.data(file="outfile",N=replicates,n=nrow(X))
	if(type=="rest") XX<-read.multi.rest.data(file="outfile",N=replicates,n=length(X))
	if(class(X)=="proseq") XX<-lapply(XX,function(x){ class(x)<-"proseq"; x })
	if(type=="sequence"||type=="morph") X<-lapply(XX,function(x,y){ rownames(x)<-y; x },y=rownames(X))
	else if(type=="rest") X<-lapply(XX,function(x,y){ names(x)<-y; x },y=names(X))
	if(!is.null(ancestors)){
		A<-readLines("outancestors")
		A<-strsplit(paste(A,collapse=""),"")[[1]]
		A<-A[A!=" "]
		m<-length(A)/replicates
		A<-lapply(1:replicates,function(i,m,x) x[1:m+(i-1)*m],m=m,x=A)
	}
	if(!is.null(mixture)){
		M<-readLines("outmixture")
		M<-strsplit(paste(M,collapse=""),"")[[1]]
		M<-M[M!=" "]
		m<-length(M)/replicates
		M<-lapply(1:replicates,function(i,m,x) x[1:m+(i-1)*m],m=m,x=M)
	}
	if(!is.null(rate.categories)){
		R<-readLines("outcategories")
		R<-strsplit(paste(R,collapse=""),"")[[1]]
		R<-R[R!=" "]
		m<-length(R)/replicates
		R<-lapply(1:replicates,function(i,m,x) as.numeric(x[1:m+(i-1)*m]),m=m,x=R)
	}
	if(!is.null(ancestors)){
		if(is.null(mixture)) X<-mapply(function(x,y) list(data=x,ancestors=y),x=X,y=A,SIMPLIFY=FALSE)
		else X<-mapply(function(x,y,z) list(data=x,ancestors=y,mixture=z),x=X,y=A,z=M,SIMPLIFY=FALSE)
	} else if(!is.null(mixture)) X<-mapply(function(x,y) list(data=x,mixture=y),x=X,y=M,SIMPLIFY=FALSE)
	if(!is.null(rate.categories)) X<-mapply(function(x,y) list(data=x,categories=y),x=X,y=R,SIMPLIFY=FALSE)
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(ancestors)) files<-c(files,"ancestors","outancestors")
		if(!is.null(rate.categories)) files<-c(files,"categories","outcategories")
		if(!is.null(mixture)) files<-c(files,"mixture","outmixture")
		cleanFiles(files)
	}
	return(X)
}

## calls clique from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rclique<-function(X,path=NULL,...){
	if(class(X)=="DNAbin") stop("you should be using Rdnacomp for DNA data.\n\n")
	if(is.vector(X)&&!is.list(X)) X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	if(is.list(X)){
		if(all(sapply(X,length)!=1)) X<-t(sapply(X,function(x) x))
		else X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	}
	if(is.data.frame(X)) X<-as.matrix(X)
	if(is.null(path)) path<-findPath("clique")
	if(is.null(path)) stop("No path provided and was not able to find path to clique")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("ancestors","infile","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r","r")
	if(hasArg(ancestral)){
		oo<-c(oo,"a")
		write(paste(ancestral,collapse=""),file="ancestors")
	} else ancestral<-NULL
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(hasArg(minimum.clique)) minimum.clique<-list(...)$minimum.clique
	else minimum.clique<-NULL
	if(!is.null(minimum.clique)) oo<-c(oo,"c",minimum.clique)	
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/clique",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(ancestral)) files<-c(files,"ancestors")
		if(!is.null(weights)) files<-c(files,"weights")
		cleanFiles(files)
	}
	return(tree)
}

## calls restml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rrestml<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("restml")
	if(is.null(path)) stop(paste("No path provided and was not able to find path to restml"))
	if(class(X)!="rest.data") stop("X should be an object of class 'rest.data'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(all.sites)) all.sites<-list(...)$all.sites
	else all.sites<-FALSE
	if(all.sites) oo<-c(oo,"a")
	if(hasArg(speedier)) speedier<-list(...)$speedier
	else speedier<-FALSE
	if((!speedier)) oo<-c(oo,"s")
	if(hasArg(global)) global<-list(...)$global
	else global<-TRUE
	if(global) oo<-c(oo,"g")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(site.length)) site.length<-list(...)$site.length
	else site.length<-6
	if(site.length!=6) oo<-c(oo,"l",site.length)
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.rest.data(X)
	system("touch outtree")
	system("touch outfile")
	temp<-system(paste(path,"/restml",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:length(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=names(X))
		cat("\n")
	}
	tree$tip.label<-names(X)[as.numeric(tree$tip.label)]
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	tree$logLik<-logLik
	return(tree)
}

## function converts string or matrix to object of class "rest.data"
## written by Liam J. Revell 2014

as.rest.data<-function(x,...){
	if(is.data.frame(x)){
		if(ncol(x)==1) x<-setNames(as.vector(t(X)),rownames(X))
		else {
			X<-as.data.frame(t(x))
			class(X)<-"rest.data"
		}	
	}
	if(is.vector(x)){
		X<-sapply(x,strsplit,split="")
		X<-as.data.frame(X)
		class(X)<-"rest.data"
	}
	if(is.matrix(x)){
		X<-as.data.frame(apply(x,1,factor))
		class(X)<-"rest.data"
	}
	attr(X,"row.names")<-NULL
	if(hasArg(nenzymes)) nenzymes<-list(...)$nenzymes
	else nenzymes<-1
	attr(X,"nenzymes")<-nenzymes
	nsites<-sapply(X,length)
	if(any(nsites!=max(nsites))) stop("All species should have the same number of sites scored.")
	attr(X,"nsites")<-min(nsites)
	return(X)
}

## S3 print method for "rest.data"
## written by Liam J. Revell 2014

print.rest.data<-function(x,printlen=6,digits=3,...){
	cat(paste(attr(x,"nsites")," restriction site scores for ",length(x)," species stored in a object of class \"rest.data\".\n\n",sep=""))
	cat(paste("All sequences of same length:",attr(x,"nsites"),"\n\n"))
	cat(paste("Number of restriction enzymes used to generate the data:",attr(x,"nenzymes"),"\n\n"))
	if(printlen<length(x))
		cat(paste("Labels:",paste(names(x)[1:min(printlen,length(x))],collapse=" "),"...\n\n"))
	else
		cat(paste("Labels:",paste(names(x)[1:min(printlen,length(x))],collapse=" "),"\n\n"))
}

## function writes rest.data to file in PHYLIP format with numbers as labels
## written by Liam J. Revell 2014

write.rest.data<-function(X,append=FALSE){
	write(paste("    ",length(X),"   ",attr(X,"nsites"),"   ",attr(X,"nenzymes"),sep=""),file="infile",append=append)
	for(i in 1:length(X)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(X[[i]],collapse=""),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
}

## calls restdist from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rrestdist<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("restdist")
	if(is.null(path)) stop("No path provided and was not able to find path to restdist")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile"))==0) return(NULL)
	oo<-c("r"); ee<-vector()
	if(hasArg(data)) data<-list(...)$data
	else data<-"sites"
	data<-tolower(data)
	if(data=="fragments") oo<-c(oo,"r")
	if(hasArg(method)) method<-list(...)$method
	else method<-"modified"
	method<-tolower(method)
	if(method=="nei/li") oo<-c(oo,"n")
	if(hasArg(gamma)){
		gamma<-list(...)$gamma
		oo<-c(oo,"g")
		ee<-c(ee,1/sqrt(gamma))
	}
	if(hasArg(kappa)){
		kappa<-list(...)$kappa
		oo<-c(oo,"t",kappa)
	}## calls gendist from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rgendist<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("gendist")
	if(is.null(path)) stop("No path provided and was not able to find path to gendist")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile"))==0) return(NULL)
	if(is.matrix(X)){
		## assumes X is a matrix of continuous character data
		N<-nrow(X)
		tips<-rownames(X)
		if(hasArg(nalleles)) nalleles<-list(...)$nalleles
		else nalleles<-rep(2,ncol(X))
		write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile")
		write(paste(nalleles,collapse=" "),file="infile",append=TRUE)
		for(i in 1:nrow(X)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			tt<-paste(sp,paste(X[i,],collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
	} else if(is.list(X)){
		## assumes X is a list of matrices containing gene frequency data
		N<-nrow(X[[1]])
		tips<-rownames(X[[1]])
		X<-lapply(X,function(x,tips) x[tips,],tips=tips)
		write(paste("    ",nrow(X[[1]]),"   ",length(X),sep=""),file="infile")
		nalleles<-sapply(X,ncol)
		write(paste(nalleles,collapse=" "),file="infile",append=TRUE)
		## verify that all rows of all X sum to 1.0
		temp<-sapply(X,rowSums)
		if(!all(round(temp,2)==1)) stop("Some of the rows of X do not sum to 1.0")
		for(i in 1:length(tips)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			dd<-vector()
			for(j in 1:length(X)) dd<-c(dd,X[[j]][i,])
			tt<-paste(sp,paste(dd,collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
	} else stop("X should be a matrix (for continuous characters) or a list (for gene frequencies)")
	oo<-c("r"); ee<-vector()
	if(hasArg(method)) method<-list(...)$method
	else method<-"nei"
	method<-tolower(method)
	if(method=="nei") oo<-c(oo,"n")
	else if(method=="cavalli-sforza") oo<-c(oo,"c")
	else if(method=="reynolds") oo<-c(oo,"r")
	else {
		cat(paste("Warning:\n  don't recognize method of type",method,".\n"))
		cat("   setting method to default type.\n\n")
		oo<-c(oo,"n")
	}
	oo<-c(oo,"y")
	system("touch outfile")
	system(paste(path,"/gendist",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
	xx<-strsplit(paste(temp,collapse=" ")," ")[[1]]
	xx<-xx[xx!=""]
	D<-matrix(NA,N,N)
	for(i in 1:N) D[i,]<-as.numeric(xx[1:N+(i-1)*(N+1)+2])
	rownames(D)<-colnames(D)<-tips
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile")
		cleanFiles(files)	
	}
	return(as.dist(D))
}

	if(hasArg(site.length)) site.length<-list(...)$site.length
	else site.length<-6
	if(site.length!=6) oo<-c(oo,"l",site.length)
	oo<-c(oo,ee,"y")
	write.rest.data(X)
	system("touch outfile")
	system(paste(path,"/restdist",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
	xx<-strsplit(paste(temp,collapse=" ")," ")[[1]]
	xx<-xx[xx!=""]
	N<-length(X)
	D<-matrix(NA,N,N)
	for(i in 1:N) D[i,]<-as.numeric(xx[1:N+(i-1)*(N+1)+2])
	rownames(D)<-colnames(D)<-names(X)
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile")
		cleanFiles(files)	
	}
	return(as.dist(D))
}

## calls kitsch from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rkitsch<-function(D,path=NULL,...){
	if(class(D)=="dist"||class(D)=="data.frame") D<-as.matrix(D)
	D<-D[rownames(D),rownames(D)]
	if(is.null(path)) path<-findPath("kitsch")
	if(is.null(path)) stop("No path provided and was not able to find path to kitsch")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(method)) method<-list(...)$method
	else method<-"fm"
	if(method=="FM"||method=="fm") method<-"fm"
	else if(method=="ME"||method=="me"){
		method<-"me"
		oo<-c(oo,"d")
	} else if(method=="LS"||method=="ls") method<-"ls"
	else {
		cat("\nWarning:\n  method not recognized - using method=\"FM\"\n")
		method="fm"
	}
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(D))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(power)){
		power<-list(...)$power
		oo<-c(oo,"p",power)
	} else if(method=="ls") oo<-c(oo,"p",0)
	if(hasArg(negative)) negative<-list(...)$negative
	else negative<-FALSE
	if(!negative) oo<-c(oo,"-")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(subreplicate)) subreplicate<-list(...)$subreplicate
	else subreplicate<-FALSE
	if(subreplicate) oo<-c(oo,"s")
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.distances(D)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/kitsch",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(D),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(D))
		cat("\n")
	}
	tree$tip.label<-rownames(D)[as.numeric(tree$tip.label)]	
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	return(tree)
}

## calls gendist from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rgendist<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("gendist")
	if(is.null(path)) stop("No path provided and was not able to find path to gendist")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile"))==0) return(NULL)
	if(is.matrix(X)){
		## assumes X is a matrix of continuous character data
		N<-nrow(X)
		tips<-rownames(X)
		if(hasArg(nalleles)) nalleles<-list(...)$nalleles
		else nalleles<-rep(2,ncol(X))
		write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile")
		write(paste(nalleles,collapse=" "),file="infile",append=TRUE)
		for(i in 1:nrow(X)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			tt<-paste(sp,paste(X[i,],collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
	} else if(is.list(X)){
		## assumes X is a list of matrices containing gene frequency data
		N<-nrow(X[[1]])
		tips<-rownames(X[[1]])
		X<-lapply(X,function(x,tips) x[tips,],tips=tips)
		write(paste("    ",nrow(X[[1]]),"   ",length(X),sep=""),file="infile")
		nalleles<-sapply(X,ncol)
		write(paste(nalleles,collapse=" "),file="infile",append=TRUE)
		## verify that all rows of all X sum to 1.0
		temp<-sapply(X,rowSums)
		if(!all(round(temp,2)==1)) stop("Some of the rows of X do not sum to 1.0")
		for(i in 1:length(tips)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			dd<-vector()
			for(j in 1:length(X)) dd<-c(dd,X[[j]][i,])
			tt<-paste(sp,paste(dd,collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
	} else stop("X should be a matrix (for continuous characters) or a list (for gene frequencies)")
	oo<-c("r"); ee<-vector()
	if(hasArg(method)) method<-list(...)$method
	else method<-"nei"
	method<-tolower(method)
	if(method=="nei") oo<-c(oo,"n")
	else if(method=="cavalli-sforza") oo<-c(oo,"c")
	else if(method=="reynolds") oo<-c(oo,"r")
	else {
		cat(paste("Warning:\n  don't recognize method of type",method,".\n"))
		cat("   setting method to default type.\n\n")
		oo<-c(oo,"n")
	}
	oo<-c(oo,"y")
	system("touch outfile")
	system(paste(path,"/gendist",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
	xx<-strsplit(paste(temp,collapse=" ")," ")[[1]]
	xx<-xx[xx!=""]
	D<-matrix(NA,N,N)
	for(i in 1:N) D[i,]<-as.numeric(xx[1:N+(i-1)*(N+1)+2])
	rownames(D)<-colnames(D)<-tips
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile")
		cleanFiles(files)	
	}
	return(as.dist(D))
}

## calls dolpenny from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rdolpenny<-function(X,path=NULL,...){
	if(is.vector(X)&&!is.list(X)) X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	if(is.list(X)){
		if(all(sapply(X,length)!=1)) X<-t(sapply(X,function(x) x))
		else X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	}
	if(is.data.frame(X)) X<-as.matrix(X)
	if(is.null(path)) path<-findPath("dolpenny")
	if(is.null(path)) stop("No path provided and was not able to find path to dolpenny")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("ancestors","infile","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(method)) method<-list(...)$method
	else method<-"dollo"
	method<-tolower(method)
	if(method=="polymorphism") oo<-c(oo,"p")
	if(hasArg(groups)){
		groups<-list(...)$groups
		oo<-c(oo,"h",groups)
	}
	if(hasArg(report)){
		report<-list(...)$report
		oo<-c(oo,"f",report)
	}
	if(hasArg(simple)) simple<-list(...)$simple
	else simple<-TRUE
	if(!simple) oo<-c(oo,"s")
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(ancestral)){
		oo<-c(oo,"a")
		ancestral<-toupper(ancestral)
		write(paste(ancestral,collapse=""),file="ancestors")
	} else ancestral<-NULL
	if(hasArg(weights)){
		weights<-list(...)$weights
		if(!any(sapply(weights,"%in%",c(0,1)))){
			cat("\n\nWarning:\n  only weights of 0 & 1 are permitted\n\n")
			weights<-NULL
		} else {
			oo<-c(oo,"w")
			write(paste(weights,collapse=""),file="weights")
		}
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dolpenny",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	xx<-strsplit(temp[ii],"  ")[[1]]
	if(class(tree)=="multiPhylo") for(i in 1:length(tree)) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
	else tree$pscore<-as.numeric(xx[length(xx)])
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(ancestral)) files<-c(files,"ancestors")
		cleanFiles(files)
	}
	return(tree)
}

## calls dollop from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2014

Rdollop<-function(X,path=NULL,...){
	if(is.vector(X)&&!is.list(X)) X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	if(is.list(X)){
		if(all(sapply(X,length)!=1)) X<-t(sapply(X,function(x) x))
		else X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	}
	if(is.data.frame(X)) X<-as.matrix(X)
	if(is.null(path)) path<-findPath("dollop")
	if(is.null(path)) stop("No path provided and was not able to find path to dollop")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("ancestors","infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(method)) method<-list(...)$method
	else method<-"dollo"
	method<-tolower(method)
	if(method=="polymorphism") oo<-c(oo,"p")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(ancestral)){
		oo<-c(oo,"a")
		ancestral<-toupper(ancestral)
		write(paste(ancestral,collapse=""),file="ancestors")
	} else ancestral<-NULL
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dollop",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	xx<-strsplit(temp[ii],"  ")[[1]]
	if(class(tree)=="multiPhylo") for(i in 1:length(tree)) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
	else tree$pscore<-as.numeric(xx[length(xx)])
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(intree)) files<-c(files,"intree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(ancestral)) files<-c(files,"ancestors")
		cleanFiles(files)
	}
	return(tree)
}

## calls penny from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rpenny<-function(X,path=NULL,...){
	if(is.vector(X)&&!is.list(X)) X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	if(is.list(X)){
		if(all(sapply(X,length)!=1)) X<-t(sapply(X,function(x) x))
		else X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	}
	if(is.data.frame(X)) X<-as.matrix(X)
	if(is.null(path)) path<-findPath("penny")
	if(is.null(path)) stop("No path provided and was not able to find path to penny")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("ancestors","infile","mixture","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(mixture)){
		oo<-c(oo,"x")
		mixture<-toupper(mixture)
		write(paste(mixture,collapse=""),file="mixture")
	} else mixture<-NULL
	if(hasArg(method)) method<-list(...)$method
	else {
		if(!is.null(mixture)) method<-"wagner"
		else method<-"mixed"
	}
	method<-tolower(method)
	if((method=="camin-sokal"||method=="wagner")&&!is.null(mixture)){
		cat("Warning: mixture provided for mixed method.\n")
		cat("         using method=\"mixed\"\n\n")
		method<-"mixed"
	}
	if(method=="camin-sokal") oo<-c(oo,"p")
	if(hasArg(groups)){
		groups<-list(...)$groups
		oo<-c(oo,"h",groups)
	}
	if(hasArg(report)){
		report<-list(...)$report
		oo<-c(oo,"f",report)
	}
	if(hasArg(simple)) simple<-list(...)$simple
	else simple<-TRUE
	if(!simple) oo<-c(oo,"s")
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(ancestral)){
		oo<-c(oo,"a")
		ancestral<-toupper(ancestral)
		write(paste(ancestral,collapse=""),file="ancestors")
	} else ancestral<-NULL
	if(hasArg(weights)){
		weights<-list(...)$weights
		if(!any(sapply(weights,"%in%",c(0,1)))){
			cat("\n\nWarning:\n  only weights of 0 & 1 are permitted\n\n")
			weights<-NULL
		} else {
			oo<-c(oo,"w")
			write(paste(weights,collapse=""),file="weights")
		}
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/penny",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	xx<-strsplit(temp[ii],"  ")[[1]]
	if(class(tree)=="multiPhylo") for(i in 1:length(tree)) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
	else tree$pscore<-as.numeric(xx[length(xx)])
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(ancestral)) files<-c(files,"ancestors")
		if(!is.null(mixture)) files<-c(files,"mixture")
		cleanFiles(files)
	}
	return(tree)
}

## calls mix from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rmix<-function(X,path=NULL,...){
	if(is.vector(X)&&!is.list(X)) X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	if(is.list(X)){
		if(all(sapply(X,length)!=1)) X<-t(sapply(X,function(x) x))
		else X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	}
	if(is.data.frame(X)) X<-as.matrix(X)
	if(is.null(path)) path<-findPath("mix")
	if(is.null(path)) stop("No path provided and was not able to find path to mix")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("ancestors","infile","intree","mixture","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(mixture)){
		oo<-c(oo,"x")
		mixture<-toupper(mixture)
		write(paste(mixture,collapse=""),file="mixture")
	} else mixture<-NULL
	if(hasArg(method)) method<-list(...)$method
	else {
		if(!is.null(mixture)) method<-"wagner"
		else method<-"mixed"
	}
	method<-tolower(method)
	if((method=="camin-sokal"||method=="wagner")&&!is.null(mixture)){
		cat("Warning: mixture provided for mixed method.\n")
		cat("         using method=\"mixed\"\n\n")
		method<-"mixed"
	}
	if(method=="camin-sokal") oo<-c(oo,"p")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(ancestral)){
		oo<-c(oo,"a")
		ancestral<-toupper(ancestral)
		write(paste(ancestral,collapse=""),file="ancestors")
	} else ancestral<-NULL
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/mix",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		if(length(ii)>1) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
		else tree$pscore<-as.numeric(xx[length(xx)])
	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(intree) files<-c(files,"intree")
		if(!is.null(ancestral)) files<-c(files,"ancestors")
		if(!is.null(mixture)) files<-c(files,"mixture")
		cleanFiles(files)
	}
	return(tree)
}

## convert phangorn phyDat to phylip.data
## written by Liam J. Revell 2013

as.phylip.data<-function(x,...){
	if(class(x)=="phyDat"&&attr(x,"type")=="USER"){
		X<-matrix(NA,length(x),length(attr(x,"index")))
		rownames(X)<-names(x)
		for(i in 1:ncol(X)){
			ii<-sapply(x,function(x,y,i) x[y[i]],y=attr(x,"index"),i=i)
			X[,i]<-attr(x,"levels")[ii]
		}
		X<-toupper(X)
		class(X)<-"phylip.data"
		return(X)
	} else {
		cat("Warning:\n  cannot convert object x to object of class 'phylip.data'.\n")
		cat("  returning NULL. Sorry!\n\n")
		return(NULL)
	}
}

## S3 print method for "phylip.data"
## written by Liam J. Revell 2013

print.phylip.data<-function(x,printlen=6,digits=3,...){
	type<-if(is.list(x)) "list" else "matrix"
	N<-if(type=="list") length(x) else nrow(x)
	cat(paste(N," character value sequences stored in a ",type,".\n\n",sep=""))
	l<-if(type=="list") sapply(x,length) else ncol(x)
	if(type=="list"){
		cat(paste("Mean sequence length:",round(mean(l),digits),"\n"))
		cat(paste("   Shortest sequence:",min(l),"\n"))
		cat(paste("    Longest sequence:",max(l),"\n\n"))
		if(N>printlen)
			cat(paste("Labels:",paste(names(x)[1:min(printlen,N)],collapse=" "),"...\n\n"))
		else
			cat(paste("Labels:",paste(names(x)[1:min(printlen,N)],collapse=" "),"\n\n"))
	} else { 
		cat(paste("All sequences of same length:",l,"\n\n"))
		if(N>printlen)
			cat(paste("Labels:",paste(rownames(x)[1:min(printlen,N)],collapse=" "),"...\n\n"))
		else
			cat(paste("Labels:",paste(rownames(x)[1:min(printlen,N)],collapse=" "),"\n\n"))
	}
	cat("Trait value composition:\n")
	ff<-summary(as.factor(x))
	print(round(ff/sum(ff),digits))
}

## calls pars from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rpars<-function(X,path=NULL,...){
	if(class(X)=="DNAbin"){ 
		cat("Warning:\n  You should be using Rdnapars for DNA data.\n\n")
		X<-as.character(X)
	}
	if(is.vector(X)&&!is.list(X)) X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	if(is.list(X)){
		if(all(sapply(X,length)!=1)) X<-t(sapply(X,function(x) x))
		else X<-t(sapply(X,function(x) strsplit(x,split="")[[1]]))
	}
	if(is.data.frame(X)) X<-as.matrix(X)
	if(is.null(path)) path<-findPath("pars")
	if(is.null(path)) stop("No path provided and was not able to find path to pars")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(thorough)) thorough<-list(...)$thorough
	else thorough<-TRUE
	if(!thorough) oo<-c(oo,"s","n")
	if(hasArg(nsave)) nsave<-list(...)$nsave
	else nsave<-100
	if(nsave!=100) oo<-c(oo,"v",nsave)
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/pars",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		if(length(ii)>1) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
		else tree$pscore<-as.numeric(xx[length(xx)])
	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	return(tree)
}

## setPath & clearPath
## written by Liam J. Revell 2013

.RphylipEnv<-new.env()
phylip.path<-NULL
setPath<-function(path) assign("phylip.path",path,envir=.RphylipEnv)
clearPath<-function() if(exists("phylip.path",envir=.RphylipEnv)) rm(phylip.path,envir=.RphylipEnv)

## calls fitch from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rfitch<-function(D,path=NULL,...){
	if(class(D)=="dist"||class(D)=="data.frame") D<-as.matrix(D)
	D<-D[rownames(D),rownames(D)]
	if(is.null(path)) path<-findPath("fitch")
	if(is.null(path)) stop("No path provided and was not able to find path to fitch")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(method)) method<-list(...)$method
	else method<-"fm"
	if(method=="FM"||method=="fm") method<-"fm"
	else if(method=="ME"||method=="me"){
		method<-"me"
		oo<-c(oo,"d")
	} else if(method=="LS"||method=="ls") method<-"ls"
	else {
		cat("\nWarning:\n  method not recognized - using method=\"FM\"\n")
		method="fm"
	}
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(D))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(power)){
		power<-list(...)$power
		oo<-c(oo,"p",power)

	} else if(method=="ls") oo<-c(oo,"p",0)
	if(hasArg(negative)) negative<-list(...)$negative
	else negative<-TRUE
	if(!negative) oo<-c(oo,"-")
	if(hasArg(global)) global<-list(...)$global
	else global<-TRUE
	if(global) oo<-c(oo,"g")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.distances(D)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/fitch",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(D),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(D))
		cat("\n")
	}
	tree$tip.label<-rownames(D)[as.numeric(tree$tip.label)]	
	if(hasArg(outgroup)){
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	return(tree)
}

## call dnainvar from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnainvar<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("dnainvar")
	if(is.null(path)) stop("No path provided and was not able to find path to dnainvar")
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(nrow(X)>4) stop("X should contain no more than 4 aligned sequences.")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2,3,4)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dnainvar",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
#	ii<-grep("total number of compatible sites is",temp)
#	for(i in 1:length(ii)){
#		xx<-strsplit(temp[ii[i]],"  ")[[1]]
#		if(length(ii)>1) tree[[i]]$compatible.sites<-as.numeric(xx[length(xx)])
#		else tree$compatible.sites<-as.numeric(xx[length(xx)])
#	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
#	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
#	else if(class(tree)=="multiPhylo"){
#		foo<-function(x,y){
#			x$tip.label<-y[as.numeric(x$tip.label)]
#			x
#		}
#		tree<-lapply(tree,foo,y=rownames(X))
#		class(tree)<-"multiPhylo"
#	}	
#	if(hasArg(outgroup)){ 
#		outgroup<-list(...)$outgroup
#		tree<-outgroup.root(tree,outgroup,quiet)
#	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile")
		if(!is.null(weights)) files<-c(files,"weights")
		cleanFiles(files)
	}
#	return(tree)
}

## call dnacomp from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnacomp<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("dnacomp")
	if(is.null(path)) stop("No path provided and was not able to find path to dnacomp")
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dnacomp",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("total number of compatible sites is",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		if(length(ii)>1) tree[[i]]$compatible.sites<-as.numeric(xx[length(xx)])
		else tree$compatible.sites<-as.numeric(xx[length(xx)])
	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	return(tree)
}

## call protdist from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rprotdist<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("protdist")
	if(is.null(path)) stop("No path provided and was not able to find path to protdist")
	if(class(X)!="proseq") stop("X should be an object of class 'proseq'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","weights","categories"))==0) return(NULL)
	oo<-c("r"); ee<-vector()
	if(hasArg(model)) model<-list(...)$model
	else model<-"JTT"
	if(model!="JTT") oo<-c(oo,rep("p",which(c("PMB","PAM","Kimura","similarity","categories")==model)))
	if(model!="Kimura"&&model!="similarity"){
		if(hasArg(gamma)){
			gamma<-list(...)$gamma
			oo<-c(oo,"g")
			ee<-c(ee,1/sqrt(gamma))
		}
	}
	if(model=="categories"){
		if(hasArg(kappa)){
			kappa<-list(...)$kappa
			oo<-c(oo,"t",kappa)
		}
		if(hasArg(bf)){
			bf<-list(...)$bf
			bf<-bf/sum(bf)
			bf<-paste(bf,collapse=" ")
			oo<-c(oo,"f",bf)
		}
		if(hasArg(genetic.code)){ 
			genetic.code<-list(...)$genetic.code
			oo<-c(oo,"u")
			genetic.code<-tolower(genetic.code)
			if(genetic.code=="universal") oo<-c(oo,"u")
			else if(genetic.code=="mitochondrial") oo<-c(oo,"m")
			else if(genetic.code=="vertebrate.mitochondrial") oo<-c(oo,"v")
			else if(genetic.code=="fly.mitochondrial") oo<-c(oo,"f")
			else if(genetic.code=="yeast.mitochondrial") oo<-c(oo,"y")
			else {
				cat(paste("Warning:\n  don't recognize genetic code of type",genetic.code,".\n"))
				cat("   setting genetic code to type 'universal'.\n\n")
				oo<-c(oo,"u")
			}
		}
		if(hasArg(categorization)){
			categorization<-list(...)$categorization
			oo<-c(oo,"a")
			categorization<-tolower(categorization)
			if(categorization=="ghb") oo<-c(oo,"g")
			else if(categorization=="chemical") oo<-c(oo,"c")
			else if(categorization=="hall") oo<-c(oo,"h")
			else {
				cat(paste("Warning:\n  don't recognize categorization of type",categorization,".\n"))
				cat("   setting categorization to default type.\n\n")
				oo<-c(oo,"g")
			}
		}
		if(hasArg(ease)){
			ease<-list(...)$ease
			oo<-c(oo,"e",ease)
		}
	}
	if(hasArg(rates)){
		rates<-list(...)$rates
		if(hasArg(rate.categories)){
			rate.categories<-list(...)$rate.categories
			write(paste(rate.categories,collapse=""),file="categories")
			ncats<-length(rates)
			rates<-paste(rates,collapse=" ")
			oo<-c(oo,"c",ncats,rates)
		} else {
			warning("cannot use rates argument without rate categories; ignoring argument rates")
			rates<-NULL
		}
	} else rates<-NULL
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")

	} else weights<-NULL
	oo<-c(oo,ee,"y")
	system("touch outfile")
	write.dna(X)
	system(paste(path,"/protdist",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
	xx<-strsplit(paste(temp,collapse=" ")," ")[[1]]
	xx<-xx[xx!=""]
	D<-matrix(NA,nrow(X),nrow(X))
	for(i in 1:nrow(X)) D[i,]<-as.numeric(xx[1:nrow(X)+(i-1)*(nrow(X)+1)+2])
	rownames(D)<-colnames(D)<-rownames(X)
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(rates)) files<-c(files,"categories")
		cleanFiles(files)
	}
	return(as.dist(D))
}

## calls protpars from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rprotpars<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("protpars")
	if(is.null(path)) stop("No path provided and was not able to find path to protpars")
	if(class(X)!="proseq") stop("X should be an object of class 'proseq'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(genetic.code)){ 
		genetic.code<-list(...)$genetic.code
		oo<-c(oo,"c")
		genetic.code<-tolower(genetic.code)
		if(genetic.code=="universal") oo<-c(oo,"u")
		else if(genetic.code=="mitochondrial") oo<-c(oo,"m")
		else if(genetic.code=="vertebrate.mitochondrial") oo<-c(oo,"v")
		else if(genetic.code=="fly.mitochondrial") oo<-c(oo,"f")
		else if(genetic.code=="yeast.mitochondrial") oo<-c(oo,"y")
		else {
			cat(paste("Warning:\n  don't recognize genetic code of type",genetic.code,".\n"))
			cat("   setting genetic code to type 'universal'.\n\n")
			oo<-c(oo,"u")
		}
	}
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/protpars",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		if(length(ii)>1) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
		else tree$pscore<-as.numeric(xx[length(xx)])
	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	return(tree)
}

## convert phangorn phyDat to proseq
## written by Liam J. Revell 2013

as.proseq<-function(x,...){
	if(class(x)=="phyDat"&&attr(x,"type")=="AA"){
		X<-matrix(NA,length(x),length(attr(x,"index")))
		rownames(X)<-names(x)
		for(i in 1:ncol(X)){
			ii<-sapply(x,function(x,y,i) x[y[i]],y=attr(x,"index"),i=i)
			X[,i]<-attr(x,"levels")[ii]
		}
		X<-toupper(X)
		class(X)<-"proseq"
		return(X)
	} else {
		cat("Warning:\n  cannot convert object x to object of class 'proseq'.\n")
		cat("  returning NULL. Sorry!\n\n")
		return(NULL)
	}
}

## read protein sequences from file
## written by Liam J. Revell 2013

read.protein<-function(file,format="fasta",...){
	X<-readLines(file)
	if(format=="fasta"){
		ii<-grep(">",X)
		nn<-gsub(">","",X[ii])
		Y<-setNames(vector("list",length=length(nn)),nn)
		ii<-cbind(ii+1,c(ii[2:length(ii)]-1,length(X)))
		for(i in 1:nrow(ii)) Y[[i]]<-strsplit(gsub(" ","",paste(X[ii[i,1]:ii[i,2]],collapse="")),"")[[1]]
		l<-sapply(Y,length)
		if(all(l==min(l))) Y<-t(sapply(Y,function(x) x))
		Y<-toupper(Y)
		class(Y)<-"proseq"
	} else if(format=="sequential"){
		xx<-strsplit(X[1]," ")[[1]]
		N<-as.numeric(xx[1])
		l<-as.numeric(xx[2])
		Y<-matrix(NA,N,l); nn<-vector()
		for(i in 1:N){
			xx<-strsplit(X[i+1]," ")[[1]]
			xx<-xx[xx!=""]
			nn[i]<-xx[1]
			Y[i,]<-strsplit(paste(xx[2:length(xx)],collapse=""),"")[[1]]
		}
		rownames(Y)<-nn
		Y<-toupper(Y)
		class(Y)<-"proseq"
	}
	return(Y)
}

## S3 print method for "proseq"
## written by Liam J. Revell 2013

print.proseq<-function(x,printlen=6,digits=3,...){
	type<-if(is.list(x)) "list" else "matrix"

	N<-if(type=="list") length(x) else nrow(x)
	cat(paste(N," protein sequences in character format stored in a ",type,".\n\n",sep=""))
	l<-if(type=="list") sapply(x,length) else ncol(x)
	if(type=="list"){
		cat(paste("Mean sequence length:",round(mean(l),digits),"\n"))
		cat(paste("   Shortest sequence:",min(l),"\n"))
		cat(paste("    Longest sequence:",max(l),"\n\n"))
		if(N>printlen)
			cat(paste("Labels:",paste(names(x)[1:min(printlen,N)],collapse=" "),"...\n\n"))
		else
			cat(paste("Labels:",paste(names(x)[1:min(printlen,N)],collapse=" "),"\n\n"))
	} else { 
		cat(paste("All sequences of same length:",l,"\n\n"))
		if(N>printlen)
			cat(paste("Labels:",paste(rownames(x)[1:min(printlen,N)],collapse=" "),"...\n\n"))
		else
			cat(paste("Labels:",paste(rownames(x)[1:min(printlen,N)],collapse=" "),"\n\n"))
	}
	cat("Amino acid composition:\n")
	ff<-summary(as.factor(x))
	print(round(ff/sum(ff),digits))
}

## calls dnamlk from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rpromlk<-function(X,path=NULL,...){
	Rproml(X,path,clock=TRUE,...)
}
	
## calls proml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rproml<-function(X,path=NULL,...){
	if(hasArg(clock)) clock<-list(...)$clock
	else clock<-FALSE
	exe<-if(clock) "promlk" else "proml"
	if(is.null(path)) path<-findPath(exe)
	if(is.null(path)) stop(paste("No path provided and was not able to find path to",exe))
	if(class(X)!="proseq") stop("X should be an object of class 'proseq'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("categories","infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r"); ee<-vector()
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(model)) model<-list(...)$model
	else model<-"JTT"
	if(model=="PMB") oo<-c(oo,"p")
	else if(model=="PAM") oo<-c(oo,rep("p",2))
	if(hasArg(rates)){
		rates<-list(...)$rates
		if(hasArg(rate.categories)){
			rate.categories<-list(...)$rate.categories
			write(paste(rate.categories,collapse=""),file="categories")
			ncats<-length(rates)
			rates<-paste(rates,collapse=" ")
			oo<-c(oo,"c",ncats,rates)
		} else {
			warning("cannot use rates argument without rate categories; ignoring argument rates")
			rates<-NULL
		}
	} else rates<-NULL
	if(hasArg(gamma)) gamma<-list(...)$gamma
	else gamma<-NULL
	if(hasArg(inv)) inv<-list(...)$inv
	else inv<-NULL
	if(hasArg(ncat)) ncat<-list(...)$ncat
	else ncat<-4
	if(!is.null(gamma)&&is.null(inv)){
		oo<-c(oo,"r")
		ee<-c(ee,1/sqrt(gamma),ncat)
	} else if(!is.null(gamma)&&!is.null(inv)){
		oo<-c(oo,"r","r")
		ee<-c(ee,1/sqrt(gamma),inv)
	}	
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(hasArg(speedier)) speedier<-list(...)$speedier
	else speedier<-FALSE
	if((!speedier)&&(!clock)) oo<-c(oo,"s")
	if(hasArg(global)) global<-list(...)$global
	else global<-TRUE
	if(global) oo<-c(oo,"g")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y",ee,"r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	temp<-system(paste(path,"/",exe,sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		if(!clock) tree<-outgroup.root(tree,outgroup,quiet)
		else cat("\nMolecular clock trees are already rooted!\n\nIgnoring argument outgroup.\n\n")
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(rates)) files<-c(files,"rates")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	tree$logLik<-logLik
	return(tree)
}

## call consense from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013, 2014

Rconsense<-function(trees,path=NULL,...){
	if(is.null(path)) path<-findPath("consense")
	if(is.null(path)) stop("No path provided and was not able to find path to consense")
	if(class(trees)!="multiPhylo") stop("trees should be an object of class 'multiPhylo'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("intree","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(method)) method<-list(...)$method
	else method<-"extended"
	if(is.numeric(method)) oo<-c(oo,rep("c",3),method)
	else if(method=="strict") oo<-c(oo,"c")
	else if(method=="majority") oo<-c(oo,rep("c",2))
	if(hasArg(outgroup)){
		outgroup<-list(...)$outgroup
		trees<-outgroup.root(trees,outgroup,quiet)
	}
	if(hasArg(rooted)) rooted<-list(...)$rooted
	else rooted<-FALSE
	if(rooted) oo<-c(oo,"r")
	if(quiet) oo<-c(oo,"2")
	oo<-c(oo,"y","r")
	tip.label<-sort(trees[[1]]$tip.label)
	trees<-lapply(trees,function(x,y){ x$tip.label<-sapply(x$tip.label,function(y,z) which(z==y),z=y,USE.NAMES=FALSE); x },y=tip.label)
	class(trees)<-"multiPhylo"
	write.tree(trees,file="intree")
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/consense",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	tree$tip.label<-tip.label[as.numeric(tree$tip.label)]
	temp<-readLines("outfile")
	if(!is.null(tree$edge.length)){
		tree$node.label<-c(NA,tree$edge.length[sapply(2:tree$Nnode+length(tree$tip.label),function(x,y) which(y==x),y=tree$edge[,2])]/length(trees))
		tree$edge.length<-NULL
	}
	if(!rooted) tree<-unroot(tree)
	if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:length(tip.label),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=tip.label)
		cat("\n")
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup) cleanFiles(c("intree","outfile","outtree"))
	return(tree)
}

## call dnapenny from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnapenny<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("dnapenny")
	if(is.null(path)) stop("No path provided and was not able to find path to dnapenny")
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(groups)){
		groups<-list(...)$groups
		oo<-c(oo,"h",groups)
	}
	if(hasArg(report)){
		report<-list(...)$report
		oo<-c(oo,"f",report)
	}
	if(hasArg(simple)) simple<-list(...)$simple
	else simple<-TRUE
	if(!simple) oo<-c(oo,"s")
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(weights)){
		weights<-list(...)$weights
		if(!any(sapply(weights,"%in%",c(0,1)))){
			cat("\n\nWarning:\n  only weights of 0 & 1 are permitted\n\n")
			weights<-NULL
		} else {
			oo<-c(oo,"w")
			write(paste(weights,collapse=""),file="weights")
		}
	} else weights<-NULL
	if(quiet) oo<-c(oo,2,3)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dnapenny",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		if(length(ii)>1) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
		else tree$pscore<-as.numeric(xx[length(xx)])
	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		cleanFiles(files)
	}
	return(tree)
}

## call dnadist from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnadist<-function(X,method=c("F84","K80","JC","LogDet"),path=NULL,...){
	method<-method[1]
	if(is.null(path)) path<-findPath("dnadist")
	if(is.null(path)) stop("No path provided and was not able to find path to dnadist")
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","weights"))==0) return(NULL)
	oo<-c("r"); ee<-vector()
	if(method!="F84") oo<-c("r",rep("d",which(c("K80","JC","LogDet","similarity")==method)))
	if(hasArg(gamma)){
		gamma<-list(...)$gamma
		oo<-c(oo,"g")
		ee<-c(ee,1/sqrt(gamma))
	}
	if(hasArg(kappa)){
		kappa<-list(...)$kappa
		oo<-c(oo,"t",kappa)
	}
	if(hasArg(rates)){
		rates<-list(...)$rates
		if(hasArg(rate.categories)){
			rate.categories<-list(...)$rate.categories
			write(paste(rate.categories,collapse=""),file="categories")
			ncats<-length(rates)
			rates<-paste(rates,collapse=" ")
			oo<-c(oo,"c",ncats,rates)
		} else {
			warning("cannot use rates argument without rate categories; ignoring argument rates")
			rates<-NULL
		}
	} else rates<-NULL
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(hasArg(bf)){
		bf<-list(...)$bf
		bf<-bf/sum(bf)
		bf<-paste(bf,collapse=" ")
		oo<-c(oo,"f",bf)
	}
	oo<-c(oo,ee,"y")


	system("touch outfile")
	write.dna(X)
	system(paste(path,"/dnadist",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
	xx<-strsplit(paste(temp,collapse=" ")," ")[[1]]
	xx<-xx[xx!=""]
	D<-matrix(NA,nrow(X),nrow(X))
	for(i in 1:nrow(X)) D[i,]<-as.numeric(xx[1:nrow(X)+(i-1)*(nrow(X)+1)+2])
	rownames(D)<-colnames(D)<-rownames(X)
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup)cleanFiles(c("infile","outfile"))
	return(as.dist(D))
}

## call treedist from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013, 2014

Rtreedist<-function(trees,method=c("branch.score","symmetric"),path=NULL,...){
	method<-method[1]
	if(is.null(path)) path<-findPath("treedist")
	if(is.null(path)) stop("No path provided and was not able to find path to treedist")
	if(hasArg(trees2)) trees2<-list(...)$trees2
	else trees2<-NULL
	N1<-if(class(trees)=="phylo") 1 else length(trees)
	if(!is.null(trees2)) N2<-if(class(trees2)=="phylo") 1 else length(trees2)
	else N2<-N1
	if(class(trees)!="multiPhylo"){
		if(!(class(trees)=="phylo"&&(class(trees2)=="phylo"||class(trees2)=="multiPhylo")))
			stop("tree should be an object of class 'phylo'")
	}
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("intree","intree2","outfile"))==0) return(NULL)
	if(class(trees)=="multiPhylo"){
		tip.label<-sort(trees[[1]]$tip.label)
		trees<-lapply(trees,function(x,y){ x$tip.label<-sapply(x$tip.label,function(y,z) which(z==y),z=y,USE.NAMES=FALSE); x },y=tip.label)
		class(trees)<-"multiPhylo"
	} else if(class(trees)=="phylo"){
		tip.label<-sort(trees$tip.label)
		trees$tip.label<-sapply(trees$tip.label,function(x,y) which(y==x),y=tip.label,USE.NAMES=FALSE)
	}
	write.tree(trees,file="intree")
	if(!is.null(trees2)){
		if(class(trees2)=="multiPhylo"){
			trees2<-lapply(trees2,function(x,y){ x$tip.label<-sapply(x$tip.label,function(y,z) which(z==y),z=y,USE.NAMES=FALSE); x },y=tip.label)
			class(trees2)<-"multiPhylo"
		} else if(class(trees2)=="phylo") trees2$tip.label<-sapply(trees2$tip.label,function(x,y) which(y==x),y=tip.label,USE.NAMES=FALSE)
		write.tree(trees2,file="intree2")
	}
	oo<-c("r")
	if(method=="symmetric") oo<-c(oo,"d")
	if(hasArg(rooted)) rooted<-list(...)$rooted
	else rooted<-FALSE
	if(rooted) oo<-c(oo,"r")
	if(quiet) oo<-c(oo,1)
	if(hasArg(distances)) distances<-list(...)$distances
	else {
		if(is.null(trees2)) distances<-"all"
		else distances<-"all.1to2"
	}
	oo<-c(oo,2)
	if(distances=="all") oo<-c(oo,"p","f")
	else if(distances=="adjacent") oo<-c(oo,"a","s")
	else if(distances=="corresponding"){
		if(is.null(trees2)){
			cat("\nWarning:")
			cat("\n  distances=\"corresponding\" not permitted for one tree object\n")
			cat("\n  defaulting to distances=\"all\"\n\n")
			distances<-"all"
			oo<-c(oo,"p","f")
		} else oo<-c(oo,"c","s")
	} else if(distances=="all.1to2"){
		if(is.null(trees2)){
			cat("\nWarning:")
			cat("\n  distances=\"all.1to2\" not permitted for one tree object\n")
			cat("\n  defaulting to distances=\"all\"\n\n")
			distances<-"all"
			oo<-c(oo,"p","f")
		} else oo<-c(oo,"l","f")
	}
	oo<-c(oo,"y")
	system("touch outfile")
	system(paste(path,"/treedist",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
	if(distances=="all"||distances=="all.1to2"){
		cc<-if(method=="symmetric") 10 else 7
		D<-matrix(NA,N1,N2)
		if(distances=="all.1to2"){
			rownames(D)<-paste(1,1:N1,sep=",")
			colnames(D)<-paste(2,1:N2,sep=",")
		} else rownames(D)<-colnames(D)<-1:N1
		nm<-ceiling(N2/cc)-1
		ii<-grep("---",temp)
		for(i in 0:nm){
			n<-min(cc,cc*(N2/cc-i))
			for(j in 1:N1){
				x<-strsplit(temp[j+ii[i+1]]," ")[[1]]
				x<-x[x!=""]
				D[j,1:n+cc*i]<-as.numeric(x[1:n+2])
			}
		}
	} else if(distances=="adjacent"||distances=="corresponding"){
		X<-read.table(file="outfile",header=FALSE,sep=" ")
		if(distances=="adjacent") D<-setNames(X[,3],paste(X[,1],X[,2],sep=","))
		else if(distances=="corresponding") D<-setNames(X[,2],X[,1])
	}
	if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:length(tip.label),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=tip.label)
		cat("\n")
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){ 
		if(is.null(trees2)) cleanFiles(c("intree","outfile"))
		else cleanFiles(c("intree","intree2","outfile"))
	}
	return(D)
}

## function to crop to first n characters a vector of strings
## written by Liam J. Revell 2013

crop<-function(x,n=1) sapply(x,function(x) strsplit(x,"")[[1]][1:n])

## function to change vector to integers
## written by Liam J. Revell 2013

to.integers<-function(x){
	types<-sort(unique(x))
	ii<-as.integer(1:length(types)-1)
	ii[sapply(x,function(x,y) which(y==x),y=types)]
}

## calls threshml from PHYLIP (Felsenstein 2013)
## written by Liam J. Revell 2013

Rthreshml<-function(tree,X,types=NULL,path=NULL,...){
	if(is.null(path)) path<-findPath("threshml")
	if(is.null(path)) stop("No path provided and was not able to find path to threshml")
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo'")
	if(is.null(types)){
		types<-sapply(X,class)
		types[types=="numeric"]<-"continuous"
		types[types%in%c("factor","character")]<-"discrete"
	}
	types<-crop(types)
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile"))==0) return(NULL)
	tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y)[1],y=rownames(X))
	## this is for the current idiosyncratic tree input file requirement of threshml
	text<-write.tree(tree)
	text<-strsplit(text,"")[[1]]

	text<-paste(paste(text[1:(length(text)-1)],collapse=""),"0.00000000;\n",sep="")
	write(text,file="intree")
	if(any(types=="c")) write.continuous(X[,types=="c"])
	if(any(types=="d")) write.dna(apply(as.matrix(X[,types=="d"]),2,to.integers),append=any(types=="c"))
	## start populating arguments
	oo<-c("r")
	if(!any(types=="d")) oo<-c(oo,"d")
	if(any(types=="c")) oo<-c(oo,"c")
	if(hasArg(burnin)){
		burnin<-list(...)$burnin
		oo<-c(oo,"b",burnin)
	}
	if(hasArg(nchain)){
		nchain<-list(...)$nchain
		oo<-c(oo,"n",nchain)
	}
	if(hasArg(ngen)){
		ngen<-list(...)$ngen
		oo<-c(oo,"s",ngen)
	}
	if(hasArg(proposal)){
		proposal<-list(...)$proposal
		oo<-c(oo,"p",proposal)
	}
	if(hasArg(lrtest)) lrtest<-list(...)$lrtest
	else lrtest<-FALSE
	if(lrtest){
		# oo<-c(oo,"t")
		cat("/nLR-test does not seem to work yet: ignoring argument lrtest.\n\n")
	}
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y",sample(seq(1,99999,by=2),1))
	system("touch outfile")
	system(paste(path,"/threshml",sep=""),input=oo,show.output.on.console=(!quiet))
	temp<-readLines("outfile")
	cc<-which(types=="c")
	dd<-which(types=="d")
	## parse covariance matrix
	ii<-grep("Covariance matrix",temp)+5
	Covariance_matrix<-matrix(NA,ncol(X),ncol(X))
	for(i in 1:ncol(X)){
		x<-strsplit(temp[i+ii]," ")[[1]]
		Covariance_matrix[i,]<-as.numeric(x[x!=""])[1:ncol(X)+1]
	}
	Covariance_matrix<-Covariance_matrix[c(cc,dd),c(cc,dd)]
	rownames(Covariance_matrix)<-colnames(Covariance_matrix)<-colnames(X)
	## parse transform matrix 1
	ii<-grep("Transform from independent variables",temp)+4
	Transform_indepvar_liab<-matrix(NA,ncol(X),ncol(X))
	for(i in 1:ncol(X)){
		x<-strsplit(temp[i+ii]," ")[[1]]
		Transform_indepvar_liab[i,]<-as.numeric(x[x!=""])[1:ncol(X)+1]
	}
	Transform_indepvar_liab<-Transform_indepvar_liab[c(cc,dd),c(cc,dd)]
	rownames(Transform_indepvar_liab)<-colnames(Transform_indepvar_liab)<-colnames(X)
	## parse variances of change
	ii<-grep("its change",temp)+1
	Var_change<-vector()
	for(i in 1:ncol(X)){
		x<-strsplit(temp[i+ii]," ")[[1]]
		Var_change[i]<-as.numeric(x[x!=""])[2]
	}
	Var_change<-Var_change[c(cc,dd)]
	names(Var_change)<-colnames(X)
	## parse transform matrix 2
	ii<-grep("Transform from liabilities or characters",temp)+4
	Transform_liab_cont<-matrix(NA,ncol(X),ncol(X))
	for(i in 1:ncol(X)){
		x<-strsplit(temp[i+ii]," ")[[1]]
		Transform_liab_cont[i,]<-as.numeric(x[x!=""])[1:ncol(X)+1]
	}
	Transform_liab_cont<-Transform_liab_cont[c(cc,dd),c(cc,dd)]
	rownames(Transform_liab_cont)<-colnames(Transform_liab_cont)<-colnames(X)
	## done parsing
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup) cleanFiles(c("infile","intree","outfile"))
	if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	return(list(Covariance_matrix=Covariance_matrix,
		Transform_indepvar_liab=Transform_indepvar_liab,
		Var_change=Var_change,
		Transform_liab_cont=Transform_liab_cont))
}

## function to write continuous characters to file
## written by Liam J. Revell 2013

write.continuous<-function(X,append=FALSE){
	if(is.vector(X)) X<-as.matrix(X)
	write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile",append=append)
	for(i in 1:nrow(X)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(X[i,],collapse=" "),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
}

## calls dnamlk from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnamlk<-function(X,path=NULL,...){
	Rdnaml(X,path,clock=TRUE,...)
}	

## clean up files
## written by Liam J. Revell 2013

cleanFiles<-function(fs){
	if(.Platform$OS.type=="windows") for(i in 1:length(fs)) system(paste("rm",fs[i],sep=" "),show.output.on.console=FALSE)
	else for(i in 1:length(fs)) system(paste("rm",fs[i],sep=" "))
}

## sets up PHYLIP in Mac OS X (based on http://evolution.gs.washington.edu/phylip/install.html)
## written by Liam J. Revell

setupOSX<-function(path=NULL){
	if(.Platform$OS.type!="unix") stop("this function is for Mac OS X only")
	if(is.null(path)){
		## check /Applications for path to PHYLIP
		ll<-list.files("/Applications/")
		ii<-grep("phylip",ll)
		if(length(ii)>0) path<-paste("/Applications/",ll[ii],sep="")
		else stop("was not able to find path to phylip installation")
	}
	if(strsplit(path,"")[length(strsplit(path,""))]=="/"){
		path<-strsplit(path,"")
		path<-paste(path[2:length(path)-1],collapse="")
	} 
	system(paste("cp ",path,"/src/linkmac ",path,"/exe/linkmac",sep=""))
	system(paste("chmod +x ",path,"/exe/linkmac",sep=""))
	system(paste("cd ",path,"/exe/\n./linkmac",sep=""))
}		

## calls neighbor from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rneighbor<-function(D,path=NULL,...){
	if(class(D)=="dist"||class(D)=="data.frame") D<-as.matrix(D)
	D<-D[rownames(D),rownames(D)]
	if(is.null(path)) path<-findPath("neighbor")
	if(is.null(path)) stop("No path provided and was not able to find path to neighbor")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(method)) method<-list(...)$method
	else method<-"nj"
	if(method=="NJ"||method=="nj") method<-"nj"
	else if(method=="UPGMA"||method=="upgma"){
		method<-"upgma"
		oo<-c(oo,"n")
	} else {
		cat("\nWarning:\n  method not recognized - using method=\"NJ\"\n")
		method="nj"
	}
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order) oo<-c(oo,"j",sample(seq(1,99999,by=2),1))
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.distances(D)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/neighbor",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(D),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(D))
		cat("\n")
	}
	tree$tip.label<-rownames(D)[as.numeric(tree$tip.label)]	
	if(hasArg(outgroup)){
		if(method=="nj"){
			outgroup<-list(...)$outgroup
			tree<-outgroup.root(tree,outgroup,quiet)
		} else { 
			cat("\nWarning:\n  outgroup rooting not permitted for method = \"upgma\"\n")
			cat("  tree already rooted\n\n")
		}
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup) cleanFiles(c("infile","outfile","outtree"))
	return(tree)
}

## write distance matrix to file in PHYLIP format
## written by Liam J. Revell 2013

write.distances<-function(D){
	write(paste("    ",nrow(D),sep=""),file="infile")
	for(i in 1:nrow(D)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
		tt<-paste(sp,paste(D[i,],collapse=" "),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
}

## attempt to find path to PHYLIP executable
## written by Liam J. Revell 2013

findPath<-function(string){
	if(exists("phylip.path",envir=.RphylipEnv)){ 
		path<-get("phylip.path",envir=.RphylipEnv)
		if(!is.null(path)) return(path)
	}
	if(.Platform$OS.type=="windows"){
		## first, check current directory
		ll<-list.files()
		ii<-grep(string,ll)
		if(length(ii)>0) 
			if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
				return(".")
		## check C:/Program Files
		ll<-list.files("C:/Program Files/")
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Program Files/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(shortPathName(dd))
		}
		## check C:/Progam Files (x86)
		ll<-list.files("C:/Program Files (x86)/")
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Program Files (x86)/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep="")))
					return(shortPathName(dd))
		}
		## check C:/Users/Username
		uu<-strsplit(getwd(),"")[[1]]
		ii<-grep("/",uu)
		uu<-paste(uu[(ii[2]+1):(ii[3]-1)],collapse="")
		ll<-list.files(paste("C:/Users/",uu,sep=""))
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Users/",uu,"/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(shortPathName(dd))
		}
		## check C:/Users/Username/Documents
		ll<-list.files(paste("C:/Users/",uu,"/Documents",sep=""))
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("C:/Users/",uu,"/Documents/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(shortPathName(dd))

		}
		return(NULL)
	} else if(.Platform$OS.type=="unix"){
		## first, check current directory
		ll<-list.files()
		ii<-grep(string,ll)
		if(length(ii)>0) 
			if(any(ll[ii]==string)) 
				return(".")
		## check /Applications
		ll<-list.files("/Applications/")
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("/Applications/",ll[ii],"/exe",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(dd)
		}
		## check /usr/local/Cellar
		ll<-list.files("/usr/local/Cellar/")
		ii<-grep("phylip",ll)
		if(length(ii)>0){
			dd<-paste("/usr/local/Cellar/",ll[ii],"/3.695/bin",sep="")
			ll<-list.files(dd)
			ii<-grep(string,ll)
			if(length(ii)>0) 
				if(any(ll[ii]==string)||any(ll[ii]==paste(string,".exe",sep=""))) 
					return(dd)
		}
		return(NULL)
	} else return(NULL)
}

## function writes DNAbin to file in PHYLIP format with numbers as labels
## written by Liam J. Revell 2013

write.dna<-function(X,append=FALSE){
	write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile",append=append)
	for(i in 1:nrow(X)){
		sp<-as.character(i)
		sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")

		tt<-paste(sp,paste(X[i,],collapse=""),collapse=" ")
		write(tt,append=TRUE,file="infile")
	}
}

## function to outgroup root
## written by Liam J. Revell 2013

outgroup.root<-function(tree,outgroup,quiet){
	if(class(tree)=="phylo") tree<-root(tree,outgroup)
	else if(class(tree)=="multiPhylo"){
		tree<-lapply(tree,root,outgroup=outgroup)
		class(tree)<-"multiPhylo"
	}
	if(!quiet){
		cat("Rooted tree(s) with the outgroup\n")
		cat("------------------------\n")
		cat(paste(paste(outgroup,collapse=", "),"\n\n"))
	}
	return(tree)
}

## function check before overwriting files
## written by Liam J. Revell 2013

file.warn<-function(gg){
	ff<-list.files()
	gg[sapply(gg,"%in%",ff)]->gg
	if(any(sapply(gg,"%in%",ff))){
		cat(paste("Warning:\n  One or more of",paste("\"",gg,"\"",sep="",collapse=", "),
		"\n  was found in your current working directory and may be overwritten\n"))
		cat("\nPress ENTER to continue or q to quit: ")
		q<-readLines(n=1)
		if(q=="q"||q=="Q") return(0) else return(1)
	} else return(1)
}

## calls dnapars from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnapars<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("dnapars")
	if(is.null(path)) stop("No path provided and was not able to find path to dnapars")
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r")
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(thorough)) thorough<-list(...)$thorough
	else thorough<-TRUE
	if(!thorough) oo<-c(oo,"s","n")
	if(hasArg(nsave)) nsave<-list(...)$nsave
	else nsave<-10000
	if(nsave!=10000) oo<-c(oo,"v",nsave)
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(hasArg(threshold)) threshold<-list(...)$threshold
	else threshold<-0
	if(threshold!=0) oo<-c(oo,"t",threshold)
	if(hasArg(transversion)) transversion<-list(...)$transversion
	else transversion<-FALSE
	if(transversion) oo<-c(oo,"n")
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y","r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	system(paste(path,"/dnapars",sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	ii<-grep("requires a total of",temp)
	for(i in 1:length(ii)){
		xx<-strsplit(temp[ii[i]],"  ")[[1]]
		if(length(ii)>1) tree[[i]]$pscore<-as.numeric(xx[length(xx)])
		else tree$pscore<-as.numeric(xx[length(xx)])
	}
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	if(class(tree)=="phylo") tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	else if(class(tree)=="multiPhylo"){
		foo<-function(x,y){
			x$tip.label<-y[as.numeric(x$tip.label)]
			x
		}
		tree<-lapply(tree,foo,y=rownames(X))
		class(tree)<-"multiPhylo"
	}	
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	return(tree)
}

## calls contml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rcontml<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("contml")
	if(is.null(path)) stop("No path provided and was not able to find path to contml")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","outfile","outtree"))==0) return(NULL)
	oo<-c("r")
	if(is.matrix(X)){
		## assumes X is a matrix of continuous character data
		write(paste("    ",nrow(X),"   ",ncol(X),sep=""),file="infile")
		for(i in 1:nrow(X)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			tt<-paste(sp,paste(X[i,],collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
		oo<-c(oo,"c")
		if(hasArg(tree)){
			oo<-c(oo,"u")
			tree<-list(...)$tree
			tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
			write.tree(tree,"intree")
			intree<-TRUE
		} else intree<-FALSE
		if(hasArg(global)) global<-list(...)$global
		else global<-TRUE
		if(global) oo<-c(oo,"g")
		if(hasArg(random.order)) random.order<-list(...)$random.order
		else random.order<-TRUE
		if(random.order){
			if(hasArg(random.addition)) random.addition<-list(...)$random.addition
			else random.addition<-10
			oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
		}
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contml",sep=""),input=oo,show.output.on.console=(!quiet))
		tree<-read.tree("outtree")
		temp<-readLines("outfile")
		logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
		temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		if(!quiet){
			cat("Translation table\n")
			cat("-----------------\n")
			temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",
				sep="")),y=rownames(X))
			cat("\n")
		}
		tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	} else if(is.list(X)){
		## assumes X is a list of matrices containing gene frequency data
		tips<-rownames(X[[1]])
		X<-lapply(X,function(x,tips) x[tips,],tips=tips)
		write(paste("    ",nrow(X[[1]]),"   ",length(X),sep=""),file="infile")
		nalleles<-sapply(X,ncol)
		write(paste(nalleles,collapse=" "),file="infile",append=TRUE)
		## verify that all rows of all X sum to 1.0
		temp<-sapply(X,rowSums)
		if(!all(round(temp,2)==1)) stop("Some of the rows of X do not sum to 1.0")
		for(i in 1:length(tips)){
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			dd<-vector()
			for(j in 1:length(X)) dd<-c(dd,X[[j]][i,])
			tt<-paste(sp,paste(dd,collapse=" "),collapse=" ")
			write(tt,append=TRUE,file="infile")
		}
		oo<-c(oo,"a")
		if(hasArg(tree)){
			oo<-c(oo,"u")
			tree<-list(...)$tree
			tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
			write.tree(tree,"intree")
			intree<-TRUE
		} else intree<-FALSE
		if(hasArg(global)) global<-list(...)$global
		else global<-TRUE
		if(global) oo<-c(oo,"g")
		if(hasArg(random.order)) random.order<-list(...)$random.order
		else random.order<-TRUE
		if(random.order){
			if(hasArg(random.addition)) random.addition<-list(...)$random.addition
			else random.addition<-10
			oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
		}
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contml",sep=""),input=oo,show.output.on.console=(!quiet))
		tree<-read.tree("outtree")
		temp<-readLines("outfile")
		logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
		temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		if(!quiet){
			cat("Translation table\n")
			cat("-----------------\n")
			temp<-lapply(1:length(tips),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",
				sep="")),y=tips)
			cat("\n")
		}
		tree$tip.label<-tips[as.numeric(tree$tip.label)]
	} else stop("X should be a matrix (for continuous characters) or a list (for gene frequencies)")
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		tree<-outgroup.root(tree,outgroup,quiet)
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)	
	}
	tree$logLik<-logLik
	return(tree)
}

## calls dnaml from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rdnaml<-function(X,path=NULL,...){
	if(hasArg(clock)) clock<-list(...)$clock
	else clock<-FALSE
	exe<-if(clock) "dnamlk" else "dnaml"
	if(is.null(path)) path<-findPath(exe)
	if(is.null(path)) stop(paste("No path provided and was not able to find path to",exe))
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("categories","infile","intree","outfile","outtree","weights"))==0) return(NULL)
	oo<-c("r"); ee<-vector()
	if(hasArg(tree)){
		oo<-c(oo,"u")
		tree<-list(...)$tree
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y),y=rownames(X))
		write.tree(tree,"intree")
		intree<-TRUE
	} else intree<-FALSE
	if(hasArg(kappa)){
		kappa<-list(...)$kappa
		oo<-c(oo,"t",kappa)
	}
	if(hasArg(bf)){
		bf<-list(...)$bf
		bf<-bf/sum(bf)
		bf<-paste(bf,collapse=" ")
		oo<-c(oo,"f",bf)
	}
	if(hasArg(rates)){
		rates<-list(...)$rates
		if(hasArg(rate.categories)){
			rate.categories<-list(...)$rate.categories
			write(paste(rate.categories,collapse=""),file="categories")
			ncats<-length(rates)
			rates<-paste(rates,collapse=" ")
			oo<-c(oo,"c",ncats,rates)
		} else {
			warning("cannot use rates argument without rate categories; ignoring argument rates")
			rates<-NULL
		}
	} else rates<-NULL
	if(hasArg(gamma)) gamma<-list(...)$gamma
	else gamma<-NULL
	if(hasArg(inv)) inv<-list(...)$inv
	else inv<-NULL
	if(hasArg(ncat)) ncat<-list(...)$ncat
	else ncat<-4
	if(!is.null(gamma)&&is.null(inv)){
		oo<-c(oo,"r")
		ee<-c(ee,1/sqrt(gamma),ncat)
	} else if(!is.null(gamma)&&!is.null(inv)){
		oo<-c(oo,"r","r")
		ee<-c(ee,1/sqrt(gamma),inv)
	}
	if(hasArg(weights)){
		oo<-c(oo,"w")
		write(paste(weights,collapse=""),file="weights")
	} else weights<-NULL
	if(hasArg(speedier)) speedier<-list(...)$speedier
	else speedier<-FALSE
	if((!speedier)&&(!clock)) oo<-c(oo,"s")
	if(hasArg(global)) global<-list(...)$global
	else global<-TRUE
	if(global) oo<-c(oo,"g")
	if(hasArg(random.order)) random.order<-list(...)$random.order
	else random.order<-TRUE
	if(random.order){
		if(hasArg(random.addition)) random.addition<-list(...)$random.addition
		else random.addition<-10
		oo<-c(oo,"j",sample(seq(1,99999,by=2),1),random.addition)
	}
	if(quiet) oo<-c(oo,2)
	oo<-c(oo,"y",ee,"r")
	write.dna(X)
	system("touch outtree")
	system("touch outfile")
	temp<-system(paste(path,"/",exe,sep=""),input=oo,show.output.on.console=(!quiet))
	tree<-read.tree("outtree")
	temp<-readLines("outfile")
	logLik<-as.numeric(strsplit(temp[grep("Ln Likelihood",temp)],"=")[[1]][2])
	temp<-lapply(temp,function(x) { cat(x); cat("\n") })
	if(!quiet){
		cat("Translation table\n")
		cat("-----------------\n")
		temp<-lapply(1:nrow(X),function(x,y) cat(paste("\t",paste(x,y[x],sep="\t"),"\n",sep="")),y=rownames(X))
		cat("\n")
	}
	tree$tip.label<-rownames(X)[as.numeric(tree$tip.label)]
	if(hasArg(outgroup)){ 
		outgroup<-list(...)$outgroup
		if(!clock) tree<-outgroup.root(tree,outgroup,quiet)
		else cat("\nMolecular clock trees are already rooted!\n\nIgnoring argument outgroup.\n\n")
	}
	if(hasArg(cleanup)) cleanup<-list(...)$cleanup
	else cleanup<-TRUE
	if(cleanup){
		files<-c("infile","outfile","outtree")
		if(!is.null(weights)) files<-c(files,"weights")
		if(!is.null(rates)) files<-c(files,"rates")
		if(intree) files<-c(files,"intree")
		cleanFiles(files)
	}
	tree$logLik<-logLik
	return(tree)
}

## function to optimize parameters of Rdnaml
## written by Liam J. Revell 2013

opt.Rdnaml<-function(X,path=NULL,...){
	if(is.null(path)) path<-findPath("dnaml")
	if(is.null(path)) stop("No path provided and was not able to find path to dnaml")
	if(class(X)!="DNAbin") stop("X should be an object of class 'DNAbin'")
	if(hasArg(tree)) tree<-list(...)$tree
	else {
		cat("\nFinding starting tree for parameter optimization\n")
		cat("------------------------------------------------\n")
		tree<-Rdnaml(X,path)
	}
	lik<-function(par,X,tree,path){
		kappa<-par[1]
		gamma<-par[2]
		bf<-par[c(3,4,5,6)]/sum(par[c(3,4,5,6)])
		ll<-Rdnaml(X,path,tree=tree,kappa=kappa,gamma=gamma,bf=bf,speedier=TRUE,random.order=FALSE,
			global=FALSE,quiet=TRUE)$logLik
		cat("\nOptimization progress\n")
		cat("-----------------------\n")
		cat(paste("kappa: ",round(kappa,5),"; gamma: ",round(gamma,5),"; bf: [",
			paste(round(bf,4),collapse=","),"]; logLik: ",round(ll,4),"\n",sep=""))
		cat("-----------------------\n")
		ll
	}
	## bounds for optimization
	if(hasArg(bounds)) bounds<-list(...)$bounds
	else bounds<-list()
	bb=list(kappa=c(0.01,20),gamma=c(0.01,20),bf=cbind(rep(0.01,4),rep(1,4)))
	bb[(namc<-names(bounds))]<-bounds
	par<-c(1,10,rep(0.25,4))
	fit<-optim(par,lik,X=X,tree=tree,path=path,method="L-BFGS-B",
		lower=c(bb$kappa[1],bb$gamma[1],bb$bf[,1]),
		upper=c(bb$kappa[2],bb$gamma[2],bb$bf[,2]),control=list(fnscale=-1))
	return(list(kappa=fit$par[1],gamma=fit$par[2],
		bf=fit$par[3:6]/sum(fit$par[3:6]),logLik=fit$value))
}

## function calls contrast from PHYLIP 3.695 (Felsenstein 2013)
## written by Liam J. Revell 2013

Rcontrast<-function(tree,X,path=NULL,...){
	if(is.null(path)) path<-findPath("contrast")
	if(is.null(path)) stop("No path provided and was not able to find path to contrast")
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo'")
	if(!is.binary.tree(tree)){
		cat("Warning:\n  Tree is not binary, resolving with branches of zero length\n")
		tree<-multi2di(tree)
	}
	if(is.vector(X)) X<-as.matrix(X)
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) if(file.warn(c("infile","intree","outfile"))==0) return(NULL)
	oo<-c("r")
	if(length(unique(rownames(X)))==nrow(X)){
		tree$tip.label<-sapply(tree$tip.label,function(x,y) which(x==y)[1],y=rownames(X))
		write.tree(tree,"intree")
		write.continuous(X)
		oo<-c(oo,"c")
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contrast",sep=""),input=oo,show.output.on.console=(!quiet))
		temp<-readLines("outfile")
		ii<-grep("Contrasts",temp)
		Contrasts<-matrix(NA,tree$Nnode,ncol(X))
		for(i in 1:tree$Nnode){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Contrasts[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Covariance",temp)
		Covariance_matrix<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Covariance_matrix[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Regressions",temp)
		Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Correlations",temp)
		Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			Correlations[i,]<-as.numeric(x[x!=""])
		}
		if(hasArg(cleanup)) cleanup<-list(...)$cleanup
		else cleanup<-TRUE
		if(cleanup) cleanFiles(c("infile","intree","outfile"))
		if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		return(list(Contrasts=Contrasts,Covariance_matrix=Covariance_matrix,
			Regressions=Regressions,Correlations=Correlations))	
	} else {
		tips<-tree$tip.label
		tree$tip.label<-1:length(tree$tip.label)
		write.tree(tree,"intree")
		write(paste("    ",length(tips),"   ",ncol(X),sep=""),file="infile")
		for(i in 1:length(tree$tip.label)){
			ii<-which(rownames(X)==tips[i])
			sp<-as.character(i)
			sp<-paste(sp,paste(rep(" ",11-nchar(sp)),collapse=""),collapse="")
			tt<-paste(sp,length(ii),collapse=" ")
			write(tt,append=TRUE,file="infile")
			for(j in 1:length(ii)) write(paste(X[ii[j],],collapse=" "),append=TRUE,file="infile")
		}
		oo<-c(oo,"w")
		if(quiet) oo<-c(oo,2)
		oo<-c(oo,"y","r")
		system("touch outfile")
		system(paste(path,"/contrast",sep=""),input=oo,show.output.on.console=(!quiet))
		temp<-readLines("outfile")
		ii<-grep("Estimate of VarA",temp)
		VarA<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarA[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Estimate of VarE",temp)[1]
		VarE<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarA Regressions",temp)
		VarA.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarA.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarA Correlations",temp)
		VarA.Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarA.Correlations[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Regressions",temp)[1]
		VarE.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Regressions",temp)[1]
		VarE.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Correlations",temp)[1]
		VarE.Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			VarE.Correlations[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Estimate of VarE",temp)[2]
		nonVa.VarE<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			nonVa.VarE[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Regressions",temp)[2]
		nonVa.VarE.Regressions<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			nonVa.VarE.Regressions[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("VarE Correlations",temp)[2]
		nonVa.VarE.Correlations<-matrix(NA,ncol(X),ncol(X))
		for(i in 1:ncol(X)){
			x<-strsplit(temp[i+ii+2]," ")[[1]]
			nonVa.VarE.Correlations[i,]<-as.numeric(x[x!=""])
		}
		ii<-grep("Log likelihood with varA",temp)
		aa<-strsplit(strsplit(temp[ii],"=")[[1]][2]," ")[[1]]
		aa<-aa[aa!=""]
		logLik<-as.numeric(sub(",","",aa[1]))
		k<-2*(ncol(X)*(ncol(X)-1)/2+ncol(X))
		ii<-grep("Log likelihood without varA",temp)
		aa<-strsplit(strsplit(temp[ii],"=")[[1]][2]," ")[[1]]
		aa<-aa[aa!=""]
		nonVa.logLik<-as.numeric(sub(",","",aa[1]))
		nonVa.k<-ncol(X)*(ncol(X)-1)/2+ncol(X)		
		ChiSq<-2*(logLik-nonVa.logLik)
		P<-pchisq(ChiSq,k-nonVa.k,lower.tail=FALSE)
		if(hasArg(cleanup)) cleanup<-list(...)$cleanup
		else cleanup<-TRUE
		if(cleanup) cleanFiles(c("infile","intree","outfile"))
		if(!quiet) temp<-lapply(temp,function(x) { cat(x); cat("\n") })
		return(list(VarA=VarA,VarE=VarE,VarA.Regressions=VarA.Regressions,
			VarA.Correlations=VarA.Correlations,
			VarE.Regressions=VarE.Regressions,
			VarE.Correlations=VarE.Correlations,
			nonVa.VarE=nonVa.VarE,
			nonVa.VarE.Regressions=nonVa.VarE.Regressions,
			nonVa.VarE.Correlations=nonVa.VarE.Correlations,
			logLik=logLik,k=k,nonVa.logLik=nonVa.logLik,
			nonVa.k=nonVa.k,P=P))
	}
}
