## =====================================================
## Translate sequences as character strings into vectors 
## (one column (variable) for each state/event)
## =====================================================

seqdecomp <- function(data, var=NULL, sep="-", miss="NA", vnames=NULL) {

	## Extracting the sequences from the data set
	seqdata <- seqxtract(data, var)

	seqdata <- as.vector(seqdata)
	nbseq <- length(seqdata)

	## Splitting the character strings
	tmp <- strsplit(seqdata, split=sep)

	## We look for the max sequence length first
	sl <- sapply(tmp,length)
	lmax <- max(sl)
	
	sdecomp <- matrix(nrow=nbseq, ncol=lmax)
	rownames(sdecomp) <- paste("[",seq(1:nbseq),"]",sep="")
	if (is.null(vnames)) 
		colnames(sdecomp) <- paste("[",seq(1:lmax),"]",sep="")
	else 
		colnames(sdecomp) <- vnames

	for (i in 1:nbseq) {
		seq <- tmp[[i]]
		seq[seq==miss] <- NA
		if (sl[i] < lmax) seq <- c(seq,rep(NA,lmax-sl[i]))
		sdecomp[i,] <- seq
		}

	return(sdecomp)

}

	
