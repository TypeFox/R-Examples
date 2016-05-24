#"as.double" <- 
#function(x, ...) UseMethod("as.double")

"as.double.gwaa.data" <- 
function(x, ...) {
	if (!is(x,"gwaa.data")) stop("data argument should be gwaa.data-class")
	x <- x@gtdata
	to <- as.double(x)
	to
}

"as.double.snp.data" <- 
function(x, ...) {
	if (!is(x,"snp.data")) stop("data argument should be snp.data-class")
	tnids <- x@nids
	tnsnps <- x@nsnps
	to <- .C("get_snps_many",as.raw(x@gtps), as.integer(tnids), as.integer(tnsnps), idata = integer(tnids*tnsnps), PACKAGE="GenABEL")$idata
	to <- replace(to,(to==0),NA)
	to <- to - 1
	dim(to) <- c(tnids,tnsnps)
	colnames(to) <- x@snpnames
	rownames(to) <- x@idnames
	to
}

"as.character.gwaa.data" <- 
function(x, ...) {
	if (!is(x,"gwaa.data")) stop("data argument should be gwaa.data-class")
	x <- x@gtdata
	to <- as.character(x)
	to
}

## old function
#"as.character.snp.data" <- 
#function(x, ...) {
#	if (!is(x,"snp.data")) stop("data argument should be snp.data-class")
#	a <- as.double(x)
#	dm <- dim(a)
#	to <- ifelse(is.na(a),NA,c("1/1","1/2","2/2")[a+1])
#	dim(to) <- dm
#	colnames(to) <- x@snpnames
#	rownames(to) <- x@idnames
#	to
#}
as.character.snp.data <- function(x,...) {
	if (!is(x,"snp.data")) stop("data argument should be snp.data-class")
	rw2ch <- alleleID.raw2char.matrix()
#
#	Very bulky fix -- now as.char.snp.coding returns true coding, used to be raw...
#
#	rect <- rw2ch[as.character(x@coding),]
#
	rect <- rw2ch[as.character(as.raw(x@coding)),]
	from <- as.double(x)
	to <- from
	if (dim(to)[2] == 1) {
		to[,1] <- ifelse(is.na(from[,1]),NA,rect[from+1])
	} else {
		for (i in 1:dim(to)[2]) {
			to[,i] <- ifelse(is.na(from[,i]),NA,rect[i,][from[,i]+1])
		}
	}
	colnames(to) <- x@snpnames
	rownames(to) <- x@idnames
	to
}

"as.genotype" <- 
function(x, ...) UseMethod("as.genotype")

"as.genotype.gwaa.data" <- 
function(x, ...) {
	if (!is(x,"gwaa.data")) stop("data argument should be gwaa.data-class")
	x <- x@gtdata
	to <- as.genotype.snp.data(x)
	to
}

"as.genotype.snp.data" <- 
function(x, ...) {
	if (!is(x,"snp.data")) stop("data argument should be snp.data-class")
	if (!require(genetics)) stop("this function requires 'genetics' package to be installed")
	gdta <- data.frame(genotype(as.character(x[,1])))
	if (x@nsnps>1) for (i in 2:x@nsnps) {
		gdta <- cbind(gdta,genotype(as.character(x[,i])))
	}
	colnames(gdta) <- x@snpnames
	rownames(gdta) <- x@idnames
#	class(gdta) <- "genotype"
	gdta
}

"as.hsgeno" <- 
function(x, ...) UseMethod("as.hsgeno")

"as.hsgeno.gwaa.data" <-
function(x, ...) {
	if (!is(x,"gwaa.data")) stop("data argument should be gwaa.data-class")
	x <- x@gtdata
	to <- as.hsgeno(x)
	to
}

"as.hsgeno.snp.data" <-
function(x, ...) {
	if (!is(x,"snp.data")) stop("data argument should be snp.data-class")
	g1 <- as.double(x[,1])
	a1 <- rep(NA,length(g1))
	a2 <- rep(NA,length(g1))
	a1 <- replace(a1,(g1==0 | g1==1),1)
	a1 <- replace(a1,(g1==2),2)
	a2 <- replace(a2,(g1==0),1)
	a2 <- replace(a2,(g1==1 | g1==2),2)
	gdta <- data.frame(a1,a2)
	if (x@nsnps>1) for (i in 2:x@nsnps) {
		g1 <- as.double(x[,i])
		a1 <- rep(NA,length(g1))
		a2 <- rep(NA,length(g1))
		a1 <- replace(a1,(g1==0 | g1==1),1)
		a1 <- replace(a1,(g1==2),2)
		a2 <- replace(a2,(g1==0),1)
		a2 <- replace(a2,(g1==1 | g1==2),2)
		gdta <- cbind(gdta,a1)
		gdta <- cbind(gdta,a2)
	}
	nams <- c()
	for (i in 1:x@nsnps) nams <- c(nams,paste(x@snpnames[i],".a1",sep=""),paste(x@snpnames[i],".a2",sep=""))
	colnames(gdta) <- nams
	rownames(gdta) <- x@idnames
#	class(gdta) <- "hsgeno"
	gdta
}

"as.raw.snp.data" <- 
function (x) {
	if (!is(x,"snp.data")) stop("data argument should be of snp.data-class")
	to <- as.raw(x@gtps)
	to
}

"as.raw.snp.mx" <-
function(x) {
	if (!is(x,"snp.mx")) stop("data argument should be of snp.mx-class")
	to <- unclass(x)
	to <- as.raw(to)
	to
}

"as.data.frame.gwaa.data" <- 
function(x, ...) {
	a <- x@phdata
	a
}

"as.character.snp.coding" <- 
function(x, ...) {
	snam <- names(x)
	rect <- alleleID.raw2char()
	x <- as.raw(x)
	out <- rect[as.character(x)]
	names(out) <- snam
	out
}

"as.character.snp.strand" <- 
function(x, ...) {
	snam <- names(x)
	x <- as.raw(x)
	tmpo <- c("u","+","-")[as.numeric(x)+1]
	names(tmpo) <- snam
	tmpo
}
