## extracts the prefix of each state in a sequence

context <- function(seqdata, L) {

	if (missing(L) || is.null(L)) {
		L <- max(seqlength(seqdata))-1
	}
	
	if (L==0) { 
		res <- matrix("e", nrow=nrow(seqdata), ncol=ncol(seqdata)) 
	} else {
		res <- matrix(nrow=nrow(seqdata), ncol=ncol(seqdata)) 
		res[,1] <- "e"
		
		for (j in ncol(seqdata):2) {
			res[,j] <- seqconc(seqdata[,max(1,j-L):(j-1)])
		}
	}

	return(res)
}



