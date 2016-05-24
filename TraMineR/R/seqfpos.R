## ===============================================================
## Search for the first occurence of a given element in a sequence
## ===============================================================

statefpos <- function (seq, state) {
	pos <- which(seq==state)
	if (length(pos)==0) fpos <- NA
	else fpos <- min(pos)
	return(fpos)
	}
	

seqfpos <- function (seqdata, state) {

	if (!inherits(seqdata,"stslist")) {
		stop("data is not a sequence object, use 'seqdef' function to create one")
	}

	fpos <- apply(seqdata,1,statefpos,state)

	return(fpos)
	}
