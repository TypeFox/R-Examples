"convert.snp.mach" <- function(pedfile,mapfile,infofile,outfile,quality=0.9, column.quality = 7, strand="+", ... ) {
	
	cat("****************************************\n")
	cat("THIS PROCDURE IS NOT RECOMMENDED FOR USE\n")
	cat("CONSIDER USE OF REGRESSION-BASED PROCEDURES\n")
	cat("FOR ANALYSIS OF IMPUTED DATA!!!\n\n")
	cat("****************************************\n")
	
	cat("Converting data to raw format...\n")
	convert.snp.ped(pedfile=pedfile,mapfile=mapfile,outfile=outfile,format="mach", wslash=T, strand=strand, ... )
	if (quality<=0) {
		cat("No quality filtering - Done.\n")
    		return(invisible(0))
	}
	
	cat("Reading info-file...\n")
	info <- read.table(file=infofile,header=TRUE)
	ncols <- dim(info)[2]
	if (ncols!=7) stop("Wrong number of columns in info-file")
	nrows <- dim(info)[1]

### read in raw data 
	cat("Loading raw data...\n")
	ifile <- file(outfile,"r")
	header <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	vver <- grep(x=header,pattern="version")
	if (length(vver)>0) {ver <- as.numeric(header[vver+1]);} else {ver <- 0;}
	if (is.na(ver)) warning("Incorrect data format version number")
	if (ver > 0) {ids <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE);}
		else {ids <- header;}
	nids <- length(ids)
	cat("   ids loaded...\n")
	mnams <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	cat("   marker names loaded...\n")
	chrom <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	chrom <- as.factor(chrom);gc(verbose=FALSE)
	cat("   chromosome data loaded...\n")
	pos <- scan(file=ifile,what=double(),nlines=1,quiet=TRUE)
	cat("   map data loaded...\n")
	if (ver==0) {
		coding <- new("snp.coding",as.raw(rep(1,length(pos))))
		strand <- new("snp.strand",as.raw(rep(0,length(pos))))
	} else {
		coding <- scan(file=ifile,what=raw(),nlines=1,quiet=TRUE)
		class(coding) <- "snp.coding"
		cat("   allele coding data loaded...\n")
		strand <- scan(file=ifile,what=raw(),nlines=1,quiet=TRUE)
		class(strand) <- "snp.strand"
		cat("   strand data loaded...\n")
	}
	nsnps <- length(mnams)
	nbytes <- ceiling(nids/4)
	rdta <- scan(file=ifile,what=raw(),quiet=TRUE)
	cat("   genotype data loaded...\n")
	close(ifile)
	dim(rdta) <- c(nbytes,nsnps)
	rdta <- new("snp.mx",rdta);gc(verbose=FALSE)

	a <- snp.data(nids=nids,rawdata=rdta,idnames=ids,snpnames=mnams,chromosome=chrom,map=pos,coding=coding,strand=strand,male=rep(0,nids))
	cat("   snp.data object created...\n")
	rm(rdta,ids,mnams,chrom,pos,coding,strand);gc(verbose=FALSE)
### end read in raw data 
	if (nrows != a@nsnps) stop("Number of SNPs in ped-file and info-file is different")
	pass <- (info[,column.quality]>=quality)
	a <- a[,pass]
	psd<-table(pass)["TRUE"]
	npsd<-table(pass)["FALSE"]
	cat("SNP filtering done --",psd,"passed quality threshold,",npsd,"did not pass...\n")
	cat("Writing data...\n")
	save.snp.data(a,genofile=outfile,human=FALSE)
	
    	return(invisible(0))
}

