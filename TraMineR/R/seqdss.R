## =======================================
## Extracts distinct states from sequences
## =======================================

seqdss <- function(seqdata, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	nbseq <- nrow(seqdata)

	sl <- seqlength(seqdata) 
	maxsl <- max(sl)

	void <- attr(seqdata, "void")
	statl <- attr(seqdata,"alphabet")
	nr <- attr(seqdata, "nr")

	trans <- matrix(void, nrow=nbseq, ncol=maxsl)

	if (with.missing) {
		statl <- c(statl, nr)
	}

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

			while (idx < sl[i] & (tmpseq[idx+1]==iseq || tmpseq[idx+1]==-99)) { 
				idx <- idx+1
			}

			## The range of the numeric alphabet 
			## obtained with seqasnum is 0..n
			if (iseq!=-99) {
				trans[i,j] <- statl[(iseq+1)]
				j <- j+1
			}
			idx <- idx+1
		}
		if (j>maxcol) {maxcol <- j}

	}
	## drop=FALSE ensures that the result is a matrix even if trans has only one row
	trans <- trans[,1:(maxcol-1), drop=FALSE]

	trans <-
          suppressMessages(
                           seqdef(trans, alphabet=alphabet(seqdata),
                                  labels=stlab(seqdata),
                                  missing=nr, right=NA,
                                  cnames=paste("ST",seq(1:(maxcol-1)),sep=""),
                                  cpal=cpal(seqdata),
                                  id=rownames(seqdata),
                                  weights=attr(seqdata, "weights")))

	return(trans)
}

