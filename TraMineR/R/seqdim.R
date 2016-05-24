## RETURNS THE DIMENSION OF ONE OR MORE SEQUENCES

seqdim <- function(seqdata) {
	## search for compressed sequences
	format <- seqfcheck(seqdata)

	if (format %in% c(":","-")) seqdata <- seqdecomp(seqdata,sep=format)

	if (is.vector(seqdata)) sdim <- c(1,length(seqdata))
	else sdim <- c(nrow(seqdata),ncol(seqdata))

	names(sdim) <- c("Rows","Columns")

	return(sdim)

	}
