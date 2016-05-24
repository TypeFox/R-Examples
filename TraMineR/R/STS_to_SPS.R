## ==============================
## Convert from STS to SPS format
## ==============================

STS_to_SPS <- function(seqdata, spsformat, 
	left=NA, right="DEL", gaps=NA, missing=NA, void="%", nr="*") {

	nbseq <- seqdim(seqdata)[1]
	maxsl <- seqdim(seqdata)[2]

	out <- matrix(NA, nrow=nbseq, ncol=maxsl)

	rownames(out) <- paste("[",seq(1:nbseq),"]",sep="")
	colnames(out) <- paste("[",seq(1:maxsl),"]",sep="")

	## Defining the format options
	prefix <- substring(spsformat$xfix,1,1)
	suffix <- substring(spsformat$xfix,2,2)
	stdursep <- spsformat$sdsep

	## PREPARING THE DATA
	seqdata <- as.matrix(seqdata)
	seqdata <- seqprep(seqdata, missing=missing, left=left, gaps=gaps, right=right, void=void, nr=nr)

	for (i in 1:nbseq) {
		idx <- 1
		j <- 1

		tmpseq <- seqdata[i,]
		sl <- TraMineR.length(tmpseq, void)
	
		while (j <= sl) {
			iseq <- tmpseq[j]

			dur <- 1
			while (j < sl & tmpseq[j+1]==iseq) { 
				dur <- dur+1
				j <- j+1
			}

			## adding suffix
			sps <- paste(prefix, iseq, stdursep, dur, suffix, sep="")

			out[i,idx] <- sps

			j <- j+1
			idx <- idx+1
		}
	}

	return(out)
}

