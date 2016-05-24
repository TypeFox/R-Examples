## ==========================================
## shifting sequence
## function used by seqformat, option to='SRS'
## ===========================================

seqshift <- function (seq, nbshift) {
	seql <- length(seq)
	
	sseq <- c(rep(NA,(seql-nbshift)), seq[1:nbshift])

	return(sseq)
}


