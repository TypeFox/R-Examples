## ==============================
## Convert from STS to DSS format
## ==============================

STS_to_DSS <- function(seqdata, 
	left=NA, right="DEL", gaps=NA, missing=NA, void="%", nr="*") {

	nbseq <- seqdim(seqdata)[1]
	maxsl <- seqdim(seqdata)[2]

	out <- matrix(NA, nrow=nbseq, ncol=maxsl)

	rownames(out) <- paste("[",seq(1:nbseq),"]",sep="")

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
			
			out[i,idx] <- iseq

			while (j < sl & tmpseq[j+1]==iseq) { 
				j <- j+1
			}
	
			j <- j+1
			idx <- idx+1
		}
	}
	return(out)
}

