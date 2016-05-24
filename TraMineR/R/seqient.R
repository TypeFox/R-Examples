## =======================
## Within Sequence Entropy
## =======================

seqient <- function(seqdata, norm=TRUE, base=exp(1), with.missing=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	statl <- attr(seqdata,"alphabet")

	if (with.missing) {
		statl <- c(statl, attr(seqdata,"nr"))
	}

	nbstat <- length(statl)

	message(" [>] computing entropy for ",nrow(seqdata)," sequences ...")

	iseqtab <- seqistatd(seqdata, with.missing=with.missing)
	
	ient <- apply(iseqtab,1,entropy, base=base)
	ient <- as.matrix(ient)
	if (norm==TRUE) {
		emax <- entropy(rep(1/nbstat,nbstat), base=base) 
		ient <- ient/emax
		}

	colnames(ient) <- "Entropy"
	rownames(ient) <- rownames(seqdata)

	return(ient)

	}	
