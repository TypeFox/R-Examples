"load.gwaa.data" <-
function(phenofile = "pheno.dat", genofile = "geno.raw",force = TRUE, makemap=FALSE, sort=TRUE, id="id") {
# check that ID and SEX are correct
	dta <- read.table(phenofile,header=TRUE,as.is=TRUE)
	coln <- names(dta)
	idColumn <- match(id,coln)
	names(dta)[idColumn] <- "id"
	if (!any(names(dta)=="id",na.rm=TRUE)) 
		stop("the filed named \"id\", containing the identifier presented in both pheno- and geno- files was not found in the phenofile")
	class(dta$id) <- "character"
	if (!any(names(dta)=="sex",na.rm=TRUE)) 
		stop("the column named \"sex\", containing the male identifier was not found in the phenofile")
#### 2.8.0!
	v <- version
	if ( as.numeric(v$major) > 2 || ((as.numeric(v$major) == 2) && (as.numeric(v$minor) >= 8.0)) ) {
		a <- table(dta$sex,useNA="ifany")
	} else {
		a <- table(dta$sex,exclude=NULL)
	}
####
	if (length(a) > 2)
		stop("column named \"sex\" contains more than 2 codes")
	if (length(a) == 1 && !(names(a)[1] == 0 || names(a)[1] == 1))
		stop("the column named \"sex\" contains 1 code which is neither 0 (=female) or 1 (=male)")
	if (length(a) == 2 && names(a)[1] != 0 && names(a)[2] != 1)
		stop("the column named \"sex\" is not coded as 0=female and 1=male")
	rm(a);gc(verbose=FALSE)
	if (any(table(dta$id)>1))
		stop("there are duplicated IDs in the phenotypic data file")
	if (any(is.na(dta$id)))
		stop("there are missing IDs in the phenotypic data file")
	if (any(is.na(dta$sex)))
		stop("there are missing sex values in the phenotypic data file")
	rownames(dta) <- dta$id
# read in genotypic data
	ifile <- file(genofile,"r")
	header <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	vver <- grep(x=header,pattern="version")
	if (length(vver)>0) {ver <- as.numeric(header[vver+1]);} else {ver <- 0;}
	if (is.na(ver)) warning("Incorrect data format version number")
	if (ver > 0) {ids <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE);}
		else {ids <- header;}
	nids <- length(ids)
	cat("ids loaded...\n")
	mnams <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	cat("marker names loaded...\n")
	chrom <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	chrom <- as.factor(chrom);gc(verbose=FALSE)
	cat("chromosome data loaded...\n")
	pos <- scan(file=ifile,what=double(),nlines=1,quiet=TRUE)
	cat("map data loaded...\n")
	if (ver==0) {
		coding <- new("snp.coding",as.raw(rep(1,length(pos))))
		strand <- new("snp.strand",as.raw(rep(0,length(pos))))
	} else {
		coding <- scan(file=ifile,what=raw(),nlines=1,quiet=TRUE)
		class(coding) <- "snp.coding"
		cat("allele coding data loaded...\n")
		strand <- scan(file=ifile,what=raw(),nlines=1,quiet=TRUE)
		class(strand) <- "snp.strand"
		cat("strand data loaded...\n")
	}
	nsnps <- length(mnams)
	nbytes <- ceiling(nids/4)
	rdta <- scan(file=ifile,what=raw(),quiet=TRUE)
	cat("genotype data loaded...\n")
	close(ifile)
	dim(rdta) <- c(nbytes,nsnps)
	rdta <- new("snp.mx",rdta);gc(verbose=FALSE)

# check errors of match between pheno and geno files
	mlst0 <- match(as.character(dta$id),ids)
	for (i in 1:length(mlst0)) {
		cid <- mlst0[i];
		if (is.na(cid)) 
			cat("person with id =",as.character(dta$id)[i],"was not found in genotypic file; excluded\n")
	}
	mlst <- match(ids,as.character(dta$id))
	oerr <- 0
	for (i in 1:length(mlst)) {
		cid <- mlst[i];
		if (is.na(cid)) {
			cat("person with id =",ids[i],"was not found in phenotypic file!!! - FATAL\n")
			oerr <- oerr + 1
		}
	}
	if (oerr) stop("fatal error. update pheno-file")
	newdta <- data.frame(dta[mlst,])
	rm(dta);gc(verbose=FALSE)

	a <- snp.data(nids=nids,rawdata=rdta,idnames=ids,snpnames=mnams,chromosome=chrom,map=pos,coding=coding,strand=strand,male=newdta$sex)
	cat("snp.data object created...\n")
	rm(rdta,ids,mnams,chrom,pos,coding,strand);gc(verbose=FALSE)

#check X chromosome markers
 	if (any(a@chromosome == "X") && any(a@male == 1) && !force) {
		xmrk <- (a@chromosome == "X")
		mlst <- (a@male == 1)
		rxm <- a[mlst,xmrk]
		Xch <- Xcheck(rxm)
		rm(rxm,xmrk,mlst);gc(verbose=FALSE)
		if (Xch$xerr) {
			cat("Wrong male X genotypes (heterozygous) found in",dim(Xch$tab)[1],"occasions\n")
			cat("Error table is saved as the output object\n")
			return(Xch$tab)
		}
	}
	if (force) cat("assignment of gwaa.data object FORCED; X-errors were not checked!\n")

# make map
	if (makemap) {
		cat("increase in map order FORCED\n")
		chun <- levels(a@chromosome)
		if (any(chun != "X")) {
			numchun <- sort(as.numeric(chun[chun!="X"]))
			gsize <- max(a@map[a@chromosome == as.character(numchun[1])])/5
			if (length(numchun)>1) {
			for (i in c(2:(length(numchun)))) {
				inc <- max(a@map[a@chromosome==as.character(numchun[i-1])]) + gsize
				a@map[a@chromosome==as.character(numchun[i])] <- a@map[a@chromosome==as.character(numchun[i])] + inc
			}
			}
			if (any(chun=="X")) {
				inc <- max(a@map[a@chromosome==as.character(numchun[length(numchun)])]) + gsize
				a@map[a@chromosome=="X"] <- a@map[a@chromosome=="X"] + inc
			}
		}
	}	

	out <- new("gwaa.data",phdata=newdta,gtdata=a)
	rm(a,newdta);gc(verbose=FALSE)
	if (sort) {
#		chr <- as.character(out@gtdata@chromosome)
#		names(chr) <- names(out@gtdata@chromosome)
#		mxC <- max(as.numeric(chr[autosomal(out@gtdata)]),na.rm=T)
#		if (any(chr=="XY")) chr <- replace(chr,(chr=="XY"),(mxC+1))
#		if (any(chr=="X")) chr <- replace(chr,(chr=="X"),(mxC+2))
#		if (any(chr=="mt")) chr <- replace(chr,(chr=="mt"),(mxC+3))
#		if (any(chr=="Y")) chr <- replace(chr,(chr=="Y"),(mxC+4))
#		chr <- as.numeric(chr)
#		ord <- order(chr,out@gtdata@map)
		ord <- sortmap.internal(out@gtdata@chromosome,out@gtdata@map)
		out <- out[,ord$ix]
	}
	out
}

