## ========================================
## Concatenates vectors of states or events
## into character strings
## ========================================

sconc <- function(seq, sep, void) {

	if (is.na(void)) 
		vi <- !is.na(seq)
	else if (!is.null(void))
		vi <- seq!=void
	else vi <- 1:length(seq)

	return(paste(seq[vi], collapse=sep))
	}

seqconc <- function (data, var=NULL, sep="-", vname="Sequence", void=NA) {

	if (inherits(data,"stslist")) {
		void <- attr(data,"void")
		cseq <- apply(data, 1, sconc, sep, void)
		cseq <- as.matrix(cseq)
		rownames(cseq) <- rownames(data)
	}
	else {
		seqdata <- seqxtract(data, var)

		if (seqdim(seqdata)[1]==1) 
			cseq <- sconc(seqdata,sep, void)
		else 
			cseq <- apply(seqdata, 1, sconc, sep, void)

		cseq <- as.matrix(cseq)

		## Rows and column names for the output
		rownames(cseq) <- paste("[",seq(1:length(cseq)),"]",sep="")
	}		
	
	colnames(cseq) <- vname

	return(cseq)
}
