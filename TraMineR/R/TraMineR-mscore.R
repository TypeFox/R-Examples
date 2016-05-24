
TraMineR.mscore <- function(seq, slength, statelist, freq) {
	mscore <- 0

	for (i in 1:slength) {
		si <- which(seq[i]==statelist)
		if (length(si)==1)
			mscore <- mscore+freq[si, i]
	}
	return(mscore)
}

