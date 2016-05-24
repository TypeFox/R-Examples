"convert.snp.text" <-
function(infile,outfile,bcast=10000) {
	ifile <- file(infile,"r")
	ofile <- file(outfile,"w")
	header <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	vver <- grep(x=header,pattern="version")
	if (length(vver)>0) {ver <- as.numeric(header[vver+1]);} else {ver <- 0;}
	if (is.na(ver)) warning("Incorrect data format version number")
	if (ver > 0) {ids <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE);}
		else {ids <- header;}
	nids <- length(ids)
	mnams <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	nmrk <- length(mnams)
	chrom <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	pos <- scan(file=ifile,what=double(),nlines=1,quiet=TRUE)
	nbytes <- ceiling(nids/4)
#	print(c("Version =",ver))
	if (ver==0) {
		coding <- as.raw(rep(1,length(pos)))
		strand <- as.raw(rep(0,length(pos)))
	} else {
		coding <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
		vec <- alleleID.char2raw()
		coding <- vec[coding]
		strand <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
		vec <- as.raw(c(0,1,2))
		names(vec) <- c("u","+","-")
		strand <- vec[strand]
	}
	cat(file=ofile,"#GenABEL raw data version 0.1\n")
	cat(file=ofile,ids,"\n")
	cat(file=ofile,mnams,"\n")
	cat(file=ofile,chrom,"\n")
	cat(file=ofile,pos,"\n")
	cat(file=ofile,coding,"\n")
	cat(file=ofile,strand,"\n")
	rdta <- raw(nbytes)
	for (i in 1:nmrk) {
		gtin <- scan(file=ifile,what=integer(),nlines=1,quiet=TRUE)
		gchk <- (gtin==0 | gtin==1 | gtin ==2 | gtin ==3)
		if (!all(gchk)) {
			cat("Wrong genotype codes:\nCODE\tID\tSNP\n")
			wlst <- which(!gchk)
			for (j in 1:length(wlst)) cat(gtin[wlst[j]],"\t",ids[wlst[j]],"\t",mnams[i],"\n")
			stop("execution terminated")
		}
		rdta <- put.snps(gtin)
		cat(file=ofile,rdta,"\n")
		if (bcast && round(i/bcast) == i/bcast) cat("Converted",i,"records...\n")
	}
	close(ifile)
	close(ofile)
}

