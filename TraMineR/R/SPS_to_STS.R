## ==============================
## Convert from SPS to STS format
## ==============================

SPS_to_STS <- function(seqdata, spsformat, nr="*") {

	nbseq <- seqdim(seqdata)[1]
	trans <- matrix("", nrow=nbseq, ncol=1)

	if (spsformat$xfix!="")
		xfix <- paste("[",spsformat$xfix,"]", sep="")
	else xfix=""
	sdsep <- spsformat$sdsep
	
	for (i in 1:nbseq) {
		tmpseq <- na.omit(seqdata[i,])

		for (s in 1:length(tmpseq)) {
			sps <- strsplit(gsub(xfix,"",tmpseq[s]), split=sdsep)[[1]]
	
			seq <- sps[1]
            if (seq==nr) seq <- NA
			dur <- as.integer(sps[2])

			if (s==1) trans[i] <- paste(trans[i],seq,sep="")
			else trans[i] <- paste(trans[i],seq,sep="-")

			if (dur>1)
				for (r in 2:dur) trans[i] <- paste(trans[i],"-",seq,sep="")
		}
	}

	trans <- seqdecomp(trans)

	return(trans)
}
