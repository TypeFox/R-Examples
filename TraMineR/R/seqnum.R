## ========================================
## Change the alphabet of a sequence object
## ========================================

seqnum <- function(seqdata, with.missing=FALSE) {
	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, see seqdef function to create one")

	alphabet.orig <- attr(seqdata,"alphabet")
	levels.orig <- levels(seqdata[,1])

	nbstat <- length(alphabet.orig)

	alphabet.new <- 0:(nbstat-1)
	levels.new <- alphabet.new

	if (with.missing) {
		## Changing missing and void codes into numerical values
		nr.new <- which(levels.orig==attr(seqdata,"nr"))-1
		levels.new <- c(levels.new, nr.new)
		attr(seqdata,"nr") <- nr.new
	}
	else 
		levels.new <- c(levels.new, attr(seqdata,"nr"))

	levels.new <- c(levels.new,attr(seqdata,"void"))

	if (length(alphabet.new)!=nbstat)
		stop("lengths of old and new alphabet are different")

	for (i in 1:seqdim(seqdata)[2])
		seqdata[,i] <- factor(seqdata[,i], levels=levels.orig, labels=levels.new)

	attr(seqdata,"alphabet") <- alphabet.new

	return(seqdata)

	}
