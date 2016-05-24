"export.impute" <- function(data,genofile="impute.gen",samplefile="impute.sample",
		strandfile="impute.strand",
		cachesizeMb=128) 
{
#	cat("beta-version of export.impute:\n\tstrand file (what is format?)\n\tsample file contains IDs of individulas\n")
	
	if (!is(data,"gwaa.data")) 
		stop("Data argumet should be of gwaa.data-class")
	if (length(levels(data@gtdata@chromosome)) > 1) 
		warning("data contains > 1 chromosome; this can not be reflected in IMPUTE format")
	if(cachesizeMb<10) {
		warning("cachesizeMB < 10, set to 10")
		cachesizeMb <- 10
	}
	# write sample file
	cat("writing sample file ...")
	miss <- (1. - perid.summary(data)$CallPP)
	gend <- data@gtdata@male
	gend[gend==0] <- 2
	samples <- data.frame(
			id=idnames(data),
			Subject_id=idnames(data),
			Missing=miss,
			Gender=gend,
			stringsAsFactors=FALSE
	)
	ord <- samples$id
	samples <- merge(samples,data@phdata,by="id",all.x=T,all.y=F)
	rownames(samples) <- samples$id
	samples <- samples[ord,]
	names(samples)[1] <- "Sample_id"
	names(samples)[4] <- "Gender"
# is this correct missing code???
	samples[is.na(samples)] <- (-9)
	write(names(samples),file=samplefile,
			ncolumns=(5+dim(data@phdata)[2]-1),append=FALSE)
	write(c(0,0,0,1,rep("P",dim(data@phdata)[2]-1)),file=samplefile,
			ncolumns=(5+dim(data@phdata)[2]-1),append=TRUE)
	write.table(samples,file=samplefile,
			col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
	rm(samples)
	gc()
	cat("... done!\n")
	
# collect info
	rsNames <- as.character(data@gtdata@snpnames)
	rsPos <- as.integer(data@gtdata@map)
	coding <- as.character(data@gtdata@coding)
	allele1=alleleID.reference()[coding]
	allele2=alleleID.effective()[coding]
	rm(coding)
	gc()
	
# write strand file
	cat("writing strand file ...\n")
	strand <- as.character(data@gtdata@strand)
	if (length(unique(strand)) > 2) warning("More then two strand types in strand file")
	tmp <- matrix(c(rsNames,rsPos,strand),ncol=3)
	write(t(tmp),file=strandfile,ncolumns=3,append=FALSE)
	rm(tmp,strand)
	cat("... done!\n")
	
# write genotype file
	cat("writing genotypes ...\n")
#	freqs <- summary(data@gtdata)[,"Q.2"]
	noutsnps <- ceiling((cachesizeMb*1024*1024)/(3*8*data@gtdata@nids))
#	print(noutsnps)
	bfrom <- seq(1,data@gtdata@nsnps,noutsnps)
	bto <- bfrom+noutsnps-1
	bto[length(bto)] <- data@gtdata@nsnps
#	print(bfrom)
#	print(bto)
	for (i in 1:length(bfrom)) {
		print(paste("writing SNPs from",bfrom[i],"to",bto[i],"..."))
		tmpdata <- data@gtdata[,c(bfrom[i]:bto[i])]
		tmp <- .Call("get_impute_snp_matrix",tmpdata@nids,tmpdata@nsnps,tmpdata@gtps) #,freqs)
		rm(tmpdata);gc()
		tmp <- cbind(rsNames[bfrom[i]:bto[i]],rsNames[bfrom[i]:bto[i]],rsPos[bfrom[i]:bto[i]],
				allele1[bfrom[i]:bto[i]],allele2[bfrom[i]:bto[i]],tmp)
# this is ugly: t(out)
		if (i == 1)
			write(t(tmp),file=genofile,ncolumns=(dim(tmp)[2]),append=FALSE)
		else
			write(t(tmp),file=genofile,ncolumns=(dim(tmp)[2]),append=TRUE)
	}
	cat("... done!\n")
}
