## RECHERCHE DE SOUS-SEQUENCES

seqpm <- function(seqdata, pattern, sep="") {

	if (!inherits(seqdata,"stslist")) {
		stop("data is not a sequence object, use 'seqdef' function to create one")
	}

	seqdata <- seqconc(seqdata,sep=sep)

	pm <- grep(pattern,seqdata)
	nbocc <- length(pm)

	message(" [>] pattern ",pattern," has been found in ", nbocc," sequences")

	res <- list(data.frame(pattern,nbocc),pm)
	names(res) <- c("MTab","MIndex")


	return(res)
	}
	
	
