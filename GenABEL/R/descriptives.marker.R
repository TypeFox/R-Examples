"descriptives.marker" <- 
function(data,snpsubset,idsubset,file,mafc,hwec,snpc,idcc,digits = 3) {
	if (is(data,"gwaa.data")) {
		if (!missing(snpsubset)) data <- data@gtdata[,snpsubset]
		if (!missing(idsubset)) data <- data@gtdata[idsubset,]
		if (missing(idsubset) & missing(snpsubset)) data <- data@gtdata
	} else if (is(data,"snp.data")) {
		if (!missing(snpsubset)) data <- data[,snpsubset]
		if (!missing(idsubset)) data <- data[idsubset,]
	} else {
		stop("data argument should be of type gwaa.data or snp.data");
	}
	if (missing(mafc)) mafc <- c(0.01,0.05,0.1,0.2)
	if (missing(hwec)) hwec <- c(0.0001,0.001,0.01,0.05)
	if (missing(snpc)) snpc <- c(.9,.95,.98,.99)
	if (missing(idcc)) idcc <- c(.9,.95,.98,.99)

	out <- list()
	s <- summary(data)
	maf <-	pmin(s$Q.2,1-s$Q.2)
	smh <- mean(2*maf*(1-maf))
	ssdh <- sd(2*maf*(1-maf))
	hwe <- s$Pexact
	cll <- s$CallRate
	rm(s);gc(verbose=FALSE)
	s <- perid.summary(data)
	rm(data);gc(verbose=FALSE)
	ill <- s[,"CallPP"]
	mhet <- mean(s[,"Het"],na.rm=T)
	sdhet <- sd(s[,"Het"],na.rm=T)
	rm(s);gc(verbose=FALSE)
	out[["Minor allele frequency distribution"]] <- catable(maf,categories=mafc,digits=digits)
	out[["Cumulative distr. of number of SNPs out of HWE, at different alpha"]] <- catable(hwe,categories=hwec,cumulative=TRUE,digits=digits)
	out[["Distribution of proportion of successful genotypes (per person)"]] <- catable(ill,categories=idcc,digits=digits)
	out[["Distribution of proportion of successful genotypes (per SNP)"]] <- catable(cll,categories=snpc,digits=digits)
	out[["Mean heterozygosity for a SNP"]] <- smh
	out[["Standard deviation of the mean heterozygosity for a SNP"]] <- ssdh
	out[["Mean heterozygosity for a person"]] <- mhet
	out[["Standard deviation of mean heterozygosity for a person"]] <- sdhet
	first <- TRUE
	if (!missing(file)) 
		for (i in names(out)) {
			if (first) {
				cat(i,"\n",file=file,append=FALSE)
				first <- FALSE
			} else {
				cat(i,"\n",file=file,append=TRUE)
			}
			if (is.matrix(out[[i]])) {
				cat("",colnames(out[[i]]),"\n",file=file,sep="\t",append=TRUE)
				write.table(out[[i]],file=file,append=TRUE,row.names=TRUE,col.names=FALSE,sep="\t")
			} else {
				cat(out[[i]],"\n",file=file,append=TRUE)
			}
		}
	out
}
