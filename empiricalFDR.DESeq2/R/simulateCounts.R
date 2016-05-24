simulateCounts <-
function(deseq.object) {
	dds=deseq.object
	sizes=sizeFactors(dds)
	counts=counts(dds)
	sums=apply(counts,1,sum)
	CVs=sqrt(dispersions(dds))
	nsamples=ncol(counts)
	message("simulating...")
	pb=txtProgressBar(0,nrow(counts))
	scounts=counts*0; deg=c()
	for (g in 1:nrow(counts)) {
		 setTxtProgressBar(pb,g)
		 if (is.na(CVs[g])){ next }
	     deg=rnorm(nsamples,sd=log(CVs[g]+1,2))
	 	 sizesDeg=sizes*(2**deg)
	 	 scounts[g,]=hist(sample(c(1:nsamples),sums[g],prob=sizesDeg,replace=TRUE),breaks=c(0:nsamples),plot=F)$counts
	}	
	simdds=DESeqDataSetFromMatrix(countData=scounts, colData=colData(dds), design=design(dds))
	simdds=estimateSizeFactors(simdds)
	return(simdds)
}
