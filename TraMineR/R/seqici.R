## ================
## Complexity index
## ================

seqici <- function(seqdata, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	## Number of transitions
	trans <- seqtransn(seqdata, with.missing=with.missing, norm=TRUE)

	## Longitudinal Entropy
	ient <- suppressMessages(
		seqient(seqdata, with.missing=with.missing, norm=TRUE)
	)

	## Complexity index
	message(" [>] computing complexity index for ",nrow(seqdata)," sequences ...")
	comp.index <- sqrt(trans * ient)

	colnames(comp.index) <- "C"

	return(comp.index)
}
