## ========================================
## Extracts states durations from sequences
## ========================================

seqdur <- function(seqdata, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, see seqdef function to create one")

	nbseq <- nrow(seqdata)
	sl <- seqlength(seqdata) 

	maxsl <- max(sl)
	trans <- matrix(nrow=nbseq, ncol=maxsl)
	rownames(trans) <- rownames(seqdata)
	colnames(trans) <- paste("DUR",1:maxsl, sep="")

	seqdatanum <- seqasnum(seqdata, with.missing=with.missing)
	if (!with.missing)
		seqdatanum[is.na(seqdatanum)] <- -99

	maxcol <- 0 
	for (i in 1:nbseq) {
		idx <- 1
		j <- 1

		tmpseq <- seqdatanum[i,]
		
		while (idx <= sl[i]) {
			iseq <- tmpseq[idx]
			dur <- 1

			while (idx < sl[i] && (tmpseq[idx+1]==iseq || tmpseq[idx+1]==-99)) { 
					idx <- idx+1
					dur <- dur+1
			}

			## The range of the numeric alphabet 
			## obtained with seqasnum is 0..n
			if (iseq!=-99) {
				trans[i,j] <- dur
				j <- j+1
			}

			idx <- idx+1
		}
		if (j>maxcol) {maxcol <- j}
	}
	## drop=FALSE ensures that the result is a matrix even if trans has only one row
	trans <- trans[,1:(maxcol-1), drop=FALSE]

	return(trans)
}

