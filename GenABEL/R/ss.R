#
#defining snp.mx -- internal class
#
#setClass("genotype",contains="data.frame",package="GenABEL")
#setClass("hsgeno",contains="data.frame",package="GenABEL")

setClass("snps.cell",contains="raw",package="GenABEL")
setClass("snp.mx",contains="matrix",package="GenABEL")
setMethod("[","snp.mx",
		function (x, i, j , drop = FALSE) {
			if (missing(j)) j=c(1:ncol(x));
			if (missing(drop)) drop = FALSE;
			if (is.logical(j)) j=which(j)
			k <- c(1:nrow(x))
			x <- x@.Data[k, j , drop = FALSE]
			if (missing(i)) i=c(1:(nrow(x)*4));
			if (is.logical(i)) i=which(i)
			x <- sset(as.raw(x),length(j),nrow(x)*4,i)
			dim(x) <- c(ceiling(length(i)/4),length(j))
			if (is.matrix(x)) 
				x <- new("snp.mx",x) 
			else 
				x <- new("snps.cell",x)
			x
		})

setMethod("show","snp.mx",
		function(object) {
			nr <- dim(object)[1]
			for (i in (1:nr)) {
				cat(as.raw(object[i,]),"\n")
			}
		})


setClass("snp.coding",contains="raw",package="GenABEL")
setMethod("[",signature(x="snp.coding",i="ANY",j="missing",drop="missing"),
		function (x,i,j,...,drop) {
			nams <- names(x)
			x <- as.raw(x)
			names(x) <- nams
			out <- new("snp.coding",x[i])
			names(out) <- nams[i]
			out
		})
setMethod("show","snp.coding",
		function (object) {
			snam <- names(object)
#		x <- alleleID.raw2char()
#		out <- x[as.character(object)]
			out <- as.raw(object)
			names(out) <- snam
			print(out)
#		cat(out,"\n");
		})


setClass("snp.strand",contains="raw",package="GenABEL")
setMethod("[",signature(x="snp.strand",i="ANY",j="missing",drop="missing"),
		function (x,i,j,...,drop) {
			nams <- names(x)
			x <- as.raw(x)
			names(x) <- nams
			out <- new("snp.strand",x[i])
			names(out) <- nams[i]
			out
		})
setMethod("show","snp.strand",
		function (object) {
			snam <- names(object)
#		tmpo <- c("u","+","-")[as.numeric(object)+1]
			tmpo <- as.raw(object)
			names(tmpo) <- snam
			print(tmpo)
#		cat(tmpo,"\n");
		})


setClass("snp.data",representation(nbytes="numeric",nids="numeric",nsnps="numeric",
				idnames="character",snpnames="character",chromosome="factor",map="numeric",
				coding="snp.coding",strand="snp.strand",
				male="numeric",gtps="ANY"),package="GenABEL")
snp.data <- function (nids,rawdata,
		idnames=as.character(c(1:nids)),
		snpnames=as.character(c(1:(length(rawdata)/ceiling(nids/4)))),
		chromosome=as.factor(rep(1,(length(rawdata)/ceiling(nids/4)))),
		map=as.double(seq(1,(length(rawdata)/ceiling(nids/4)))),
		coding=as.raw(rep(1,length(rawdata)/ceiling(nids/4))),
		strand=as.raw(rep(0,length(rawdata)/ceiling(nids/4))),
		male=rep(0,nids)
) {
	nbytes <- ceiling(nids/4)
	nsnps <- length(rawdata)/nbytes
	if (nsnps != round(nsnps)) stop("wrong number of ids")
	names(map) <- snpnames
	names(chromosome) <- snpnames
	names(coding) <- snpnames
	names(strand) <- snpnames
	names(male) <- idnames
	t <- as.raw(rawdata)
	dim(t) <- c(nbytes,nsnps)
	t <- new("snp.mx",t)
	a <- new("snp.data",nbytes=nbytes,nids=nids,nsnps=nsnps,gtps=t,
			idnames=idnames,snpnames=snpnames,
			chromosome=chromosome,map=map,
			coding=coding,
			strand=strand,
			male=male)
	a
}

setMethod("[","snp.data",
		function(x, i, j, drop) {
			if (missing(j)) j=c(1:x@nsnps);
			if (missing(i)) i=c(1:x@nids);
			if (missing(drop)) drop = FALSE;
			if (is.logical(i)) i=which(i)
			if (is.logical(j)) j=which(j)
			if (is.character(i)) {
				tmp <- i
				i=match(i,x@idnames)
				if (any(is.na(i))) stop(paste("following IDs were not found:",tmp[which(is.na(i))],"\n"))
			}
			if (is.character(j)) {
				tmp <- j
				j=match(j,x@snpnames)
				if (any(is.na(j))) stop(paste("following SNPs were not found:",tmp[which(is.na(j))],"\n"))
			}
			if (length(i) > x@nids  || max(i) > x@nids ) stop("i out of range")
			if (length(j) > x@nsnps || max(j) > x@nsnps) stop("j out of range")
			if (length(i) <= 0) stop("i out of range (== 0)")
			if (length(j) <= 0) stop("j out of range (== 0)")
			a <- new("snp.data")
			a@nids <- length(i)
			a@nsnps <- length(j)
			a@nbytes <- ceiling(a@nids/4.)
			a@snpnames <- x@snpnames[j]
			a@map <- x@map[j]
			a@chromosome <- as.factor(as.character(x@chromosome[j]))
			names(a@chromosome) <- a@snpnames
			a@coding <- new("snp.coding",x@coding[j])
			names(a@coding) <- a@snpnames
			a@strand <- new("snp.strand",x@strand[j])
			names(a@strand) <- a@snpnames
			a@idnames <- x@idnames[i]
			a@male <- x@male[i]
			names(a@male) <- a@idnames
			a@gtps <- x@gtps[i,j]
			a
		})
#setMethod("summary","snp.data",
#	function(object,...) {
#)

setMethod("show","snp.data",
		function(object) {
			cat("@nids =",object@nids,"\n")
			cat("@nsnps =",object@nsnps,"\n")
			cat("@nbytes =",object@nbytes,"\n")
			cat("@idnames =",object@idnames,"\n")
			cat("@snpnames =",object@snpnames,"\n")
			cat("@chromosome =",object@chromosome,"\n")
			cat("@coding = ",as.raw(object@coding),"\n")
			cat("@strand = ",as.raw(object@strand),"\n")
			cat("@map =",object@map,"\n")
			cat("@male =",object@male,"\n")
			cat("@gtps = \n");
			show(object@gtps)
		})

setMethod(
		f = "dim",
		signature = "snp.data",
		definition = function(x)
		{
			return(c(nids(x),nsnps(x)))
		}
);

setMethod(
		f = "dimnames",
		signature = "snp.data",
		definition = function(x)
		{
			return(list(idnames(x),snpnames(x)))
		}
);


setClass("gwaa.data",representation(phdata="data.frame",gtdata="snp.data"),package="GenABEL")

setMethod("show","gwaa.data",
		function (object) {
			show(object@phdata)
			show(object@gtdata)
		})

setMethod("[","gwaa.data",
		function( x, i, j, drop) {
			if (missing(j)) j=c(1:x@gtdata@nsnps);
			if (missing(i)) i=c(1:x@gtdata@nids);
			if (missing(drop)) drop = FALSE;
			if (is.logical(i)) i=which(i)
			if (is.logical(j)) j=which(j)
			if (is.character(i)) {
				tmp <- i
				i=match(i,x@gtdata@idnames)
				if (any(is.na(i))) stop(paste("following IDs were not found:",tmp[which(is.na(i))],"\n"))
			}
			if (is.character(j)) {
				tmp <- j
				j=match(j,x@gtdata@snpnames)
				if (any(is.na(j))) stop(paste("following SNPs were not found:",tmp[which(is.na(j))],"\n"))
			}
			if (length(i) > x@gtdata@nids  || max(i) > x@gtdata@nids ) stop("i out of range")
			if (length(j) > x@gtdata@nsnps || max(j) > x@gtdata@nsnps) stop("j out of range")
			if (length(i) <= 0) stop("i out of range (== 0)")
			if (length(j) <= 0) stop("j out of range (== 0)")
			a <- x@gtdata[i,j]
			b <- x@phdata[i,]
			out <- new("gwaa.data",phdata=b,gtdata=a)
			out
		})


#setClass("scan.gwaa",contains="list",package="GenABEL")
setClass("scan.gwaa.2D",contains="list",package="GenABEL")








#
# scan.gwaa class
# (c) 2010 Yurii Aulchenko, EMCR
#
#

setClass(
		Class = "scan.gwaa",
		representation(
				results = "data.frame",
				dimresults = "integer",
				annotation = "data.frame",
				lambda = "list",
				idnames = "character",
				call = "call",
				family = "character"
		),
		package = "GenABEL"
);


setMethod(
		f = "initialize",
		signature = "scan.gwaa",
		definition = function(.Object,results,
				annotation = NULL, lambda = NULL, idnames = NULL, 
				call = NULL, family = NULL)
		{
			.Object@results <- results
			.Object@dimresults <- dim(results)
			if (!is.null(annotation)) {.Object@annotation <- annotation}
			else {.Object@annotation <- NULL}
			if (!is.null(lambda)) {.Object@lambda <- lambda}
			else {.Object@lambda <- NULL}
			if (!is.null(idnames)) {.Object@idnames <- idnames}
			else {.Object@idnames <- NULL}
			if (!is.null(call)) {.Object@call <- call}
			else {.Object@call <- NULL}
			if (!is.null(family)) {.Object@family <- family}
			else {.Object@family <- NULL}
			return(.Object)
		}
);

# replace standard methods

setMethod(
		f = "show",
		signature = "scan.gwaa",
		definition = function(object)
		{
			cat("***** 'scan.gwaa' object *****\n")
			if (!is.null(object@call)) {cat("*** Produced with:\n");print(getcall(object));}
			if (!is.null(object@family)) cat("*** Test used:",getfamily(object),"\n");
			if (!is.null(object@idnames)) cat("*** no. IDs used:",length(object@idnames),"(",
						idnames(object)[1:min(3,length(object@idnames))],", ... )\n")
			if (!is.null(object@lambda)) {
				cat("*** Lambda:",lambda(object)$est,"\n");
				#	cat("*** Lambda:\n");
				#	print(lambda(object));
			}
			cat("*** Results table contains",object@dimresults[1],"rows and",object@dimresults[2],"columns\n")
			if (object@dimresults[1]>10) {
				res <- object@results[1:10,]
				cat("*** Output for 10 first rows is:\n")
				print(res)
				cat("   ...\n")
			} else {
				res <- results(object)
				print(res)
			}
			cat("___ Use 'results(object)' to get complete results table ___\n")
		}
);

setMethod(
		f = "[",
		signature = "scan.gwaa",
		definition = function(x,i,j,drop)
		{
			res <- results(x)
			if (missing(j)) drop <- FALSE
			return(res[i,j,drop=drop])
		}
);

setMethod(
		f = "dim",
		signature = "scan.gwaa",
		definition = function(x)
		{
			return(x@dimresults)
		}
);

setMethod(
		f = "dimnames",
		signature = "scan.gwaa",
		definition = function(x)
		{
			return(dimnames(results(x)))
		}
);

#
# set new generics
#

setGeneric('is.scan.gwaa', function(x) standardGeneric('is.scan.gwaa'))

setMethod('is.scan.gwaa', signature(x='scan.gwaa'),
		function(x) return(TRUE))

# accessing slots of 'scan.gwaa'

setGeneric(
		name = "results",
		def = function(object, showGCP = TRUE ) {standardGeneric("results");}
);

setMethod(
		f = "results",
		signature = "scan.gwaa",
		definition = function(object, showGCP = TRUE ) 
		{
			if (!is.null(object@annotation)) {
				tmp <- cbind(object@annotation,object@results)
				rownames(tmp) <- rownames(object@results)
			} else {
				tmp <- object@results
			}
			if (!is.null(lambda(object)) & !any(names(tmp)=="Pc1df")) {
				chi2.1df_corr <- tmp$chi2.1df/lambda(object)$estimate
				tmp$Pc1df <- pchisq(chi2.1df_corr, df=1, lower.tail=FALSE)
			}
			return(tmp)
		}
);

setGeneric(
		name = "annotation",
		def = function(object) {standardGeneric("annotation");}
);
setMethod(
		f = "annotation",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(object@annotation)
		}
);
setMethod(
		f = "annotation",
		signature = "snp.data",
		definition = function(object) 
		{
			cod <- coding(object)
			A1 <- substr(cod,1,1)
			A2 <- substr(cod,2,2)
#			res <- data.frame(Chromosome=chromosome(object),Position=map(object),
#					Strand=strand(object),A1=A1,A2=A2,stringsAsFactors = FALSE)
# should save space
			res <- data.frame(Chromosome=chromosome(object),Position=map(object),
					Strand=strand(object),A1=A1,A2=A2,stringsAsFactors = TRUE)
			rownames(res) <- snpnames(object)
			res
		}
);
setMethod(
		f = "annotation",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(annotation(object@gtdata))
		}
);



setGeneric(
		name = "getfamily",
		def = function(object) {standardGeneric("getfamily");}
);

setMethod(
		f = "getfamily",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(object@family)
		}
);


setGeneric(
		name = "getcall",
		def = function(object) {standardGeneric("getcall");}
);

setMethod(
		f = "getcall",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(object@call)
		}
);


setGeneric(
		name = "lambda",
		def = function(object) {standardGeneric("lambda");}
);

setMethod(
		f = "lambda",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(object@lambda)
		}
);


setGeneric(
		name = "idnames",
		def = function(object) {standardGeneric("idnames");}
);
setMethod(
		f = "idnames",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(idnames(object@gtdata))
		}
);
setMethod(
		f = "idnames",
		signature = "snp.data",
		definition = function(object) 
		{
			return(object@idnames)
		}
);
setMethod(
		f = "idnames",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(object@idnames)
		}
);


setGeneric(
		name = "snpnames",
		def = function(object) {standardGeneric("snpnames");}
);
setMethod(
		f = "snpnames",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(snpnames(object@gtdata))
		}
);
setMethod(
		f = "snpnames",
		signature = "snp.data",
		definition = function(object) 
		{
			return(object@snpnames)
		}
);
setMethod(
		f = "snpnames",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(rownames(object@results))
		}
);

setGeneric(
		name = "nids",
		def = function(object) {standardGeneric("nids");}
);
setMethod(
		f = "nids",
		signature = "snp.data",
		definition = function(object) 
		{
			return(object@nids)
		}
);
setMethod(
		f = "nids",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(nids(object@gtdata))
		}
);
setMethod(
		f = "nids",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(length(idnames(object)))
		}
);

setGeneric(
		name = "nsnps",
		def = function(object) {standardGeneric("nsnps");}
);
setMethod(
		f = "nsnps",
		signature = "snp.data",
		definition = function(object) 
		{
			return(object@nsnps)
		}
);
setMethod(
		f = "nsnps",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(nsnps(object@gtdata))
		}
);
setMethod(
		f = "nsnps",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(dim(results(object))[1])
		}
);

setGeneric(
		name = "map",
		def = function(object) {standardGeneric("map");}
);
setMethod(
		f = "map",
		signature = "snp.data",
		definition = function(object) 
		{
			return(object@map)
		}
);
setMethod(
		f = "map",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(map(object@gtdata))
		}
);
setMethod(
		f = "map",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			mp <- annotation(object)[,"Position"]
			names(mp) <- snpnames(object)
			return(mp)
		}
);

setGeneric(
		name = "chromosome",
		def = function(object) {standardGeneric("chromosome");}
);
setMethod(
		f = "chromosome",
		signature = "snp.data",
		definition = function(object) 
		{
			chr <- as.character(object@chromosome)
			names(chr) <- names(object@chromosome)
			return(chr)
		}
);
setMethod(
		f = "chromosome",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(chromosome(object@gtdata))
		}
);
setMethod(
		f = "chromosome",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			chr <- as.character(annotation(object)[,"Chromosome"])
			names(chr) <- snpnames(object) 
			return(chr)
		}
);


setGeneric(
		name = "strand",
		def = function(object) {standardGeneric("strand");}
);
setMethod(
		f = "strand",
		signature = "snp.data",
		definition = function(object) 
		{
			return(as.character(object@strand))
		}
);
setMethod(
		f = "strand",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(strand(object@gtdata))
		}
);
setMethod(
		f = "strand",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(annotation(object)[,"Strand"])
		}
);

setGeneric(
		name = "strand<-",
		def = function(x,value) {standardGeneric("strand<-");}
);
setMethod(
		f = "strand<-",
		signature = "snp.data",
		definition = function(x,value) 
		{
			x <- patch_strand(data=x,snpid=snpnames(x),strand=value)
			return(x)
		}
);
setMethod(
		f = "strand<-",
		signature = "gwaa.data",
		definition = function(x,value) 
		{
			strand(x@gtdata) <- value
			return(x)
		}
);



setGeneric(
		name = "coding",
		def = function(object) {standardGeneric("coding");}
);
setMethod(
		f = "coding",
		signature = "snp.data",
		definition = function(object) 
		{
			return(as.character(object@coding))
		}
);
setMethod(
		f = "coding",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(coding(object@gtdata))
		}
);
setMethod(
		f = "coding",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			tmp <- annotation(object)[,c("A1","A2")]
			return(paste(tmp[,1],tmp[,2],sep=""))
		}
);


setGeneric(
		name = "coding<-",
		def = function(x,value) {standardGeneric("coding<-");}
);
setMethod(
		f = "coding<-",
		signature = "snp.data",
		definition = function(x,value) 
		{
			rawVal <- as.raw(alleleID.char2raw()[value])
			names(rawVal) <- x@snpnames
			x@coding <- new("snp.coding",rawVal)
			if (any(is.na(coding(x)))) {
				cat("wrong coding value, should be one of ")
				cat(names(alleleID.char2raw()),"\n")
				stop()
			}
			return(x)
		}
);
setMethod(
		f = "coding<-",
		signature = "gwaa.data",
		definition = function(x,value) 
		{
			coding(x@gtdata) <- value
			return(x)
		}
);


setGeneric(
		name = "refallele",
		def = function(object) {standardGeneric("refallele");}
);
setMethod(
		f = "refallele",
		signature = "snp.data",
		definition = function(object) 
		{
			return(substr(as.character(object@coding),1,1))
		}
);
setMethod(
		f = "refallele",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(substr(coding(object@gtdata),1,1))
		}
);
setMethod(
		f = "refallele",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(annotation(object)[,"A1"])
		}
);


setGeneric(
		name = "effallele",
		def = function(object) {standardGeneric("effallele");}
);
setMethod(
		f = "effallele",
		signature = "snp.data",
		definition = function(object) 
		{
			return(substr(as.character(object@coding),2,2))
		}
);
setMethod(
		f = "effallele",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(substr(coding(object@gtdata),2,2))
		}
);
setMethod(
		f = "effallele",
		signature = "scan.gwaa",
		definition = function(object) 
		{
			return(annotation(object)[,"A2"])
		}
);



setGeneric(
		name = "male",
		def = function(object) {standardGeneric("male");}
);
setMethod(
		f = "male",
		signature = "snp.data",
		definition = function(object) 
		{
			return(object@male)
		}
);
setMethod(
		f = "male",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(male(object@gtdata))
		}
);


setGeneric(
		name = "gtdata",
		def = function(object) {standardGeneric("gtdata");}
);
setMethod(
		f = "gtdata",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(object@gtdata)
		}
);


setGeneric(
		name = "phdata",
		def = function(object) {standardGeneric("phdata");}
);
setMethod(
		f = "phdata",
		signature = "gwaa.data",
		definition = function(object) 
		{
			return(object@phdata)
		}
);


setGeneric(
		name = "phdata<-",
		def = function(x,value) {standardGeneric("phdata<-");}
);
setMethod(
		f = "phdata<-",
		signature = "gwaa.data",
		definition = function(x,value) 
		{
			x@phdata <- value
			return(x)
		}
);

