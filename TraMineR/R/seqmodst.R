## ====================
## Modal state sequence
## ====================

seqmodst <- function(seqdata, weighted=TRUE, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, see seqdef function to create one", call.=FALSE)

	slength <- ncol(seqdata)
	statl <- alphabet(seqdata)
	cnames <- colnames(seqdata)

	if (with.missing) {
		statl <- c(statl, attr(seqdata,"nr"))
	}

	## State distribution
	freq <- seqstatd(seqdata, weighted, with.missing)$Frequencies

	ctype <- matrix(nrow=1, ncol=slength)
	stfreq <- matrix(nrow=1, ncol=slength)
	colnames(stfreq) <- cnames
	rownames(stfreq) <- "Freq."
	colnames(ctype) <- cnames

	## Constructing the transversal modal sequence
	for (i in 1:slength) {
		smax <- which(freq[,i]==max(freq[,i]))[1]
		stfreq[,i] <- freq[smax,i]
		ctype[,i] <- statl[smax]
	}
	
	res <- suppressMessages(seqdef(ctype, alphabet=alphabet(seqdata),
		missing=attr(seqdata,"nr"), nr=attr(seqdata,"nr"),
		left=NA, gaps=NA, right=NA,
		labels=stlab(seqdata),
		cpal=cpal(seqdata), missing.color=attr(seqdata,"missing.color"),
		xtstep=attr(seqdata, "xtstep")))

	nbocc <- length(seqfind(res, seqdata))

	## Distance to modal state sequence
	## if (dist)
	##	dist.modst <- seqdist(seqdata, refseq=res, ...)
	## else
	##	dist.modst <- NULL

	class(res) <- c("stslist.modst", class(res))

	attr(res, "Frequencies") <- stfreq
	## attr(res, "Distances") <- dist.modst
	attr(res, "Occurrences") <- nbocc

	## Weights
	weights <- attr(seqdata, "weights")

	if (!weighted || is.null(weights)) {
		weights <- rep(1.0, nrow(seqdata))
	}
	if (all(weights==1))
		weighted <- FALSE

	attr(res, "nbseq") <- sum(weights)
	attr(res,"weighted") <- weighted

	return(res)
 }
