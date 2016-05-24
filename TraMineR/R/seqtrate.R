## =========================
## Computes transition rates
## =========================

seqtrate <- function(seqdata, statl=NULL, time.varying=FALSE, weighted=TRUE, lag=1, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist")) { 
		stop(" [!] seqdata is NOT a sequence object, see seqdef function to create one")
	}

	## State list if not specified
	if (is.null(statl)){
		statl <- attr(seqdata,"alphabet")
	}
	if(with.missing){
		statl <- c(statl, attr(seqdata,"nr"))
	}
	nr <- attr(seqdata,"nr")
	void <- attr(seqdata,"void")

	## Weights
	weights <- attr(seqdata, "weights")
	
	if (!weighted || is.null(weights)) {
		weights <- rep(1, nrow(seqdata))
	}
	nbetat <- length(statl)

	sdur <- ncol(seqdata)

	seqdata <- as.matrix(seqdata)
	
	if (lag<0) {
		alltransition <- (abs(lag)+1):sdur
	} else {
		alltransition <- 1:(sdur-lag)
	}
	numtransition <- length(alltransition)
	
	## =================================
	## time varying transition rates
	## =================================
	if (time.varying) {
		message(" [>] computing time varying transition rates for states ",paste(statl,collapse="/")," ...",sep="")
		## Dimension names
		dnames <- list()
		dnames[[1]] <- paste("[",statl," ->]",sep="")
		dnames[[2]] <- paste("[-> ",statl,"]",sep="")
		dnames[[3]] <- colnames(seqdata)[alltransition]
		
		tmat <- array(0, dim=c(nbetat, nbetat, numtransition))
		dimnames(tmat) <- dnames
		for (sl in alltransition) {
			missingcond <- seqdata[,sl+lag]!=void
			if(!with.missing){
				missingcond <- missingcond & seqdata[,sl+lag]!=nr
			}
			for (x in 1:nbetat) {
				colxcond <- seqdata[,sl]==statl[x]
				PA <- sum(weights[colxcond & missingcond])
					
				if (PA == 0) {
					tmat[x,,sl] <- 0
				}
				else {
					for (y in 1:nbetat) {
						PAB <- sum(weights[colxcond & seqdata[,sl+lag]==statl[y]])
						tmat[x,y,sl] <- PAB/PA
					}
				}
			}
		}
	}
	## =================================
	## Non time varying transition rates
	## =================================
	else {
		message(" [>] computing transition rates for states ",paste(statl,collapse="/")," ...",sep="")
		tmat <- matrix(0, nrow=nbetat, ncol=nbetat)
		row.names(tmat) <- paste("[",statl," ->]",sep="")
		colnames(tmat) <- paste("[-> ",statl,"]",sep="")
		missingcond <- seqdata[,alltransition+lag]!=void
		if(!with.missing){
			missingcond <- missingcond & seqdata[,alltransition+lag]!=nr
		}
		for (x in 1:nbetat) {
			## Count
			PA <- 0
			colxcond <- seqdata[, alltransition]== statl[x] 
			if(numtransition>1){
				PA <- sum(weights * rowSums(colxcond & missingcond))
			}
			else{
				PA <- sum(weights * (colxcond & missingcond))
			}
			if (PA==0){
				 tmat[x, ] <- 0
			}
			else {
				for (y in 1:nbetat) {
					if(numtransition>1){
						PAB <- sum(weights * rowSums(colxcond & seqdata[,alltransition+lag]==statl[y]))
					} else{
						PAB <- sum(weights * (colxcond & seqdata[,alltransition+lag]==statl[y]))
					}
					tmat[x,y] <- PAB/PA
				}
			}
		}
	}
	return(tmat)
}


