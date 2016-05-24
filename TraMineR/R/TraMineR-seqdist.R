## ====================================================
## Computing distances to a reference sequence
## Internal function called by seqdist
## ====================================================




TraMineR.seqdist.refseq <- function(seqdata, method, refseq,
	norm, indel, sm, alphsize, nd, dseq, slength, mcorr,with.missing){
	
	## Getting refseq
	##User specified
	if (inherits(refseq,"stslist") && nrow(refseq)==1) {
		compseq <- refseq
		message(" [>] using (external) sequence ",
        	suppressMessages(seqformat(compseq, from="STS", to="SPS", compressed=TRUE)), " as reference")
	}
	## Most frequent sequence as reference
	else if (refseq==0) {
		mfseq <- seqtab(seqdata, tlim=1)
		message(" [>] using the most frequent sequence as reference: ",
			suppressMessages(seqformat(mfseq, from="STS", to="SPS", compressed=TRUE)))
		idxmfseq <- suppressMessages(seqfind(mfseq, seqdata))
		message(" [>] the most frequent sequence appears ", length(idxmfseq), " times")
		compseq <- seqdata[idxmfseq[1],]
	}
	## Indice of sequence given as reference
	else if (is.numeric(refseq) & refseq>0) {
		compseq <- seqdata[refseq,]
		message(" [>] using sequence ", refseq,": ",
   			suppressMessages(seqformat(compseq, from="STS", to="SPS", compressed=TRUE)), " as reference")
	} else {
		stop("[!] invalid reference sequence", call.=FALSE)
	}
	## Length of compseq
	lcompseq <- seqlength(compseq)
	## Vector of distance
	m <- vector(mode="numeric", length=nd)
	compseq <- seqasnum(seqnum(compseq), with.missing=with.missing)
	if (method=="OM") {
		for (i in 1:nd) {
			m[i] <- levenshtein(dseq[i,], slength[i], compseq, lcompseq, indel,sm,alphsize,norm)
		}
	} else if (method=="LCP") {
		for (i in 1:nd) {
			m[i] <- LCPdist(dseq[i,], slength[i], compseq, lcompseq, norm)
		}
	} else if (method=="RLCP") {
		## reverse order of compseq ignoring missing values at the end
		compseq <- rev(compseq[1:lcompseq])
		for (i in 1:nd) {
			m[i] <- LCPdist(rev(dseq[i,1:slength[i]]),slength[i],compseq,lcompseq,norm)
		}
	} else if (method=="LCS") {
		for (i in 1:nd) {
			m[i] <- LCSdist(dseq[i,],slength[i],compseq,lcompseq,norm)
		}
	} else if (method == "DHD" || method == "HAM") {
		for (i in 1:nd) {
			m[i] <- DHDdist(dseq[i,], slength[i], compseq, lcompseq, indel, sm, norm)
		}
	}
	## Constructing the final distance vector
	mcorr <- match(seqconc(seqdata),seqconc(dseq))
	distances <- m[mcorr]
	names(distances) <- NULL
	return(distances)
}


TraMineR.seqdist.all <- function(seqdata, method,
	norm, indel, sm, alphsize, nd, dseq, slength, mcorr, optimized){
	
	magicSeq <- order(mcorr)
	magicIndex <- c(unique(rank(mcorr, ties.method="min")), nrow(seqdata)+1)-1

	if (method=="OM") {
		## One for OM, 2 for LCP
   		disttype <- as.integer(1)
		if (optimized) {
			disttype <- as.integer(0)
		}
	}
	else if (method=="LCP") {
		disttype <- as.integer(2) ## One for OM, 2 for LCP
		sm <- 0
		indel <- 0
	}
	else if (method=="RLCP") {
		disttype <- as.integer(3) ## One for OM, 2 for LCP
		sm <- 0
		indel <- 0
	}
	else if (method=="DHD") {
		disttype <- as.integer(4) ## 4 for DHD or HAM
	}
 	
	distances <- .Call(TMR_cstringdistance,
		as.integer(dseq),
		as.integer(dim(dseq)),
		as.integer(slength),
		as.double(indel),
		as.integer(alphsize),
		as.double(sm),
		as.integer(norm),
		as.integer(magicIndex),
		as.integer(magicSeq),
		disttype)

	## Setting some attributes for the dist object
 	class(distances) <- "dist"
	attr(distances,"Size") <- length(magicSeq)
 	attr(distances,"method") <- method
 	attr(distances, "Labels") <- dimnames(seqdata)[[1]]
 	attr(distances, "Diag") <- FALSE
 	attr(distances, "Upper") <- FALSE
	return(distances)
}
