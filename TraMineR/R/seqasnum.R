## TRANSFORMATION DE SEQUENCES DE CARACTERES EN SEQUENCES NUMERIQUES

seqasnum <- function(seqdata, with.missing=FALSE) {

	mnum <- matrix(NA,nrow=seqdim(seqdata)[1],ncol=seqdim(seqdata)[2])

	rownames(mnum) <- rownames(seqdata)
	colnames(mnum) <- colnames(seqdata)

	statl <- attr(seqdata,"alphabet")

	if (with.missing)
		statl <- c(statl, attr(seqdata, "nr"))

	for (i in 1:length(statl)) {
		mnum[seqdata==statl[i]] <- i-1
		}

	return(mnum)

	}
