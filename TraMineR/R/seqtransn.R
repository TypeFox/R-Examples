## =====================================
## Number of transitions in the sequence
## =====================================

seqtransn <- function(seqdata, with.missing=FALSE, norm=FALSE, pweight=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	## Number of transitions
	dss <- seqdss(seqdata, with.missing=with.missing)
	dssl <- seqlength(dss)
	nbseq <- nrow(dss)

	if (pweight) {
		tr <- seqtrate(seqdata)
		dss.num <- seqasnum(dss)+1
		trans <- matrix(0, nrow=nbseq, ncol=1)
		rownames(trans) <- rownames(seqdata)

		for (i in 1:nbseq) {
			if (dssl[i]>1) {
				for (j in 2:dssl[i]) {
					trans[i] <- trans[i] + (1-tr[dss.num[i,j-1], dss.num[i,j]])
				}
			}
		}
	}
	else {
		trans <- dssl-1
		if (any(dssl==0)) {
			trans[dssl==0] <- 0
		}
	}

	if (norm) {
		seql <- seqlength(seqdata)
		trans <- trans/(seql-1)
		if (any(seql==1)) {
			trans[seql==1] <- 0
		}
	}

	colnames(trans) <- "Trans."

	return(trans)
}

trans.pweight <- function(seqdata, tr) {
	res <- 0

	for (i in 2:seqlength(seqdata)) {
		res <- res + (1-tr[seqdata[i-1], seqdata[i]])
	}

	return(res)
}
