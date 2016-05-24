## ====================================================
## Computing distances between sequences
## Available metrics (method):
## OM = optimal matching
## LCP = Longest Common Prefix (Elzinga)
## LCPnorm = Longest Common Prefix normalized (Elzinga)
## LCS = Longest Common Subsequence (Elzinga)
## LCSold = Long Common Subsequence (R code) (Elzinga)
## ====================================================

seqdistold <- function(seqdata, method, refseq=NULL, 
	norm=FALSE, indel=1, sm=NA,
	with.miss=FALSE) {

	## ======
	## CHECKS
	## ======
	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, use 'seqdef' function to create one")

	metlist <- c("OM","LCP", "LCS")
	if (!method %in% metlist) 
		stop("Method must be one of: ", paste(metlist,collapse=" "))

	n <- seqdim(seqdata)[1]
	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)
	message(" [>] ",n," sequences with ", alphsize,
		" distinct events/states (", paste(alphabet,collapse="/"),")")

	if (with.miss) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
		alphsize <- length(alphabet)
		message(" [>] including missing value as additional state" )
	}
	else 
		if (any(seqdata==attr(seqdata,"nr")))
			stop("found missing values in sequences, please set 'with.miss' option to nevertheless compute distances")

	## Checking if substitution cost matrix contains values for each state
	if (method=="OM") 
		if (nrow(sm)!=alphsize | ncol(sm)!=alphsize)
			stop("size of substitution cost matrix must be ",alphsize,"x", alphsize)

	## Reference sequence
	if (!is.null(refseq)) {
		if (refseq==0) {
			mfseq <- row.names(seqtab(seqdata,tlim=1))
			message(" [>] using most frequent sequence as reference: ", mfseq)

			idxmfseq <- suppressMessages(which(seqformat(seqdata, 
				from='STS', to='SPS', compressed=TRUE)==mfseq))
			message(" [>] most frequent sequence appears ", length(idxmfseq), " times")

			compseq <- seqdata[idxmfseq[1],]

		}
		if (refseq>0) {
			compseq <- seqdata[refseq,]
			message(" [>] using sequence ",refseq, seqformat(compseq,from="STS",to="SPS")," as reference")
		}
		lcompseq <- sum(!is.na(compseq))
		distmat <- FALSE
	} 
	else
		distmat <- TRUE

	## ==============
	## Preparing data
	## ==============
	seqdata <- seqnum(seqdata, with.missing=with.miss)

	## Selecting distinct sequences only and saving the indexes 
	dseq <- unique(seqdata)
	mcorr <- match(seqconc(seqdata),seqconc(dseq))

	## nd <- seqdim(dseq)[1]
	nd <- nrow(dseq)
	message(" [>] ", nd," distinct sequences")

	if (distmat==FALSE) {
		m <- vector(mode="numeric", length=nd)
		compseq <- seqasnum(seqnum(compseq))
	}
	else {
		m <- matrix(nrow=nd, ncol=nd)
		diag(m) <- 0
	}

	## l <- ncol(dseq)
	slength <- seqlength(dseq)

	dseq <- seqasnum(dseq, with.missing=with.miss)

	message(" [>] min/max sequence length: ",min(slength),"/",max(slength))

	debut <- Sys.time()

	message(" [>] computing distances using ",method, appendLF=FALSE)
	if (norm==TRUE) 
		message(" normalized metric ...",appendLF =FALSE) 
	else 
		message(" metric ...",appendLF =FALSE)
	
	## Function and arguments 
	if (method=="OM") {
		if (distmat==FALSE) {
			for (i in 1:nd) 
				m[i] <- levenshtein(dseq[i,], slength[i], compseq, lcompseq, indel,sm,alphsize,norm)
		} else {
			for(i in 1:(nd-1)) {
				l1 <- slength[i]
				seq1 <- dseq[i,1:l1]
				for(j in 1:(nd-i)) {
					l2 <- slength[i+j]
					seq2 <- dseq[i+j,1:l2]	
					m[i+j,i] <- levenshtein(seq1,l1,seq2,l2,indel,sm,alphsize,norm)
				}
			}
		}
	}

	if (method=="LCP") {
		if (distmat==FALSE) {
			for (i in 1:nd) 
				m[i] <- LCPdist(dseq[i,],slength[i],compseq,lcompseq,norm)
		} else {
			for(i in 1:(nd-1)) {
				l1 <- slength[i]
				seq1 <- dseq[i,1:l1]
				for(j in 1:(nd-i)) {
					l2 <- slength[i+j]
					seq2 <- dseq[i+j,1:l2]
					m[i+j,i] <- LCPdist(seq1,l1,seq2,l2,norm)
				}
			}
		}
	}

     if (method=="LCS") {
		if (distmat==FALSE) {
			for (i in 1:nd) 
				m[i] <- LCSdist(dseq[i,],slength[i],compseq,lcompseq,norm)
		} else {
			for(i in 1:(nd-1)) {
				l1 <- slength[i]
				seq1 <- dseq[i,1:l1]
				for(j in 1:(nd-i)) {
					l2 <- slength[i+j]
					seq2 <- dseq[i+j,1:l2]
					m[i+j,i] <- LCSdist(seq1,l1,seq2,l2,norm)
				}
			}
		}	
	}

	fin <- Sys.time()
	message(" (",round(difftime(fin,debut,units="mins"),2)," minutes)")

	## =========================================================
	## Constructing the final distance matrix with all sequences
	## =========================================================
	if (distmat==TRUE) {
		message(" [>] creating distance matrix ... ",appendLF =FALSE)
		debut2 <- Sys.time()

		m[upper.tri(m)] <- t(m)[upper.tri(m)]
		md <- matrix(nrow=n, ncol=n)
		md <- m[mcorr,mcorr]

		dimnames(md) <- list(rownames(seqdata),rownames(seqdata))
				
		fin <- Sys.time()
		message("(",round(difftime(fin,debut2,units="mins"),2)," minutes)")
	} 
	else {
		mcorr <- match(seqconc(seqdata),seqconc(dseq))
		md <- m[mcorr]
		names(md) <- NULL
	}
		
	return(md)
}

