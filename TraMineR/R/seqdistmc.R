## multichannel distances

seqdistmc <- function(channels, method, norm=FALSE, indel=1, sm=NULL,
	with.missing=FALSE, full.matrix=TRUE, link="sum", cval=2, miss.cost=2, cweight=NULL ) {
	
	## Checking arguments
	nchannels <- length(channels)
	if (nchannels < 2) {
		stop("[!] please specify at least two channels")
	}
	if (is.null(cweight)) {
		cweight <- rep(1, nchannels)
	}
	numseq <- sapply(channels,nrow)
	if(any(numseq!=numseq[1])) {
		stop(" [!] sequence objects have different numbers of rows")
	}
	numseq <- numseq[1]
	message(" [>] ", nchannels, " channels with ", numseq, " sequences")
	## Actually LCP and RLCP are not included
	metlist <- c("OM", "LCS", "DHD", "HAM")
	if (!method %in% metlist) {
		stop(" [!] method must be one of: ", paste(metlist, collapse=" "), call.=FALSE)
	}
	## We handle LCS as OM
	if (method=="LCS") {
		method <- "OM"
		sm <- "CONSTANT"
		indel <- 1
		cval <- 2
		miss.cost <- 2
	}
	timeVarying <- method %in% c("DHD")
	## Indels and sm only apply for OM
	## Correct number of arguments
	## Generating default substitution arguments for DHD and HAM
	if (is.null(sm)) {
		costmethod <- "CONSTANT"
		if (method == "DHD") {
			costmethod <- "TRATE"
		}
		sm <- rep(costmethod, nchannels)
	}
	else if (length(sm)==1 && sm %in% c("CONSTANT", "TRATE")){
		sm <- rep(sm, nchannels)
	}
	if (length(indel)==1) {
		indel <- rep(indel, nchannels)
	}
	## Checking correct numbers of info per channel
	if ((length(indel)!= nchannels) ||
		(length(sm)!= nchannels) ||
		(length(cweight)!= nchannels)) {
		stop(" [!] you should supply one weight, one substitution matrix and one indel per channel")
	}
	## indels
	indel_list <- numeric(length=nchannels)
	## subsitution matrix
	substmat_list <- list()
	## alphabet for each channel
	alphabet_list <- list()
	## alphabet size per channel
	alphsize_list <-list()
	## seqlenth of each channels
	maxlength_list <- numeric(length=nchannels)
	
	## ============================================================
	## Building and checking substitution matrix per channel
	## ============================================================
	for (i in 1:nchannels) {
		## Sequence object
		if (!inherits(channels[[i]],"stslist")) {
			stop(" [!] channel ", i, " is not a state sequence object, use 'seqdef' function to create one", call.=FALSE)
		}
		alphabet_list[[i]] <- attr(channels[[i]],"alphabet")
		## Checking missing values
		if (with.missing) {
			alphabet_list[[i]] <- c(alphabet_list[[i]],attr(channels[[i]],"nr"))
			message(" [>] including missing value as an additional state" )
		}
		else {
			if (any(channels[[i]]==attr(channels[[i]],"nr"))) {
				stop(" [!] found missing values in channel ", i, ", please set 'with.missing=T' to nevertheless compute distances")
			}
		}
		alphsize_list[[i]] <- length(alphabet_list[[i]])
		## Storing number of columns
		maxlength_list[i] <- ncol(channels[[i]])
		indel_list[i] <- indel[i]
		## Substitution matrix generation method is given
		if	(is.character(sm[[i]])) {
			message(" [>] computing substitution cost matrix for channel ", i)
			substmat_list[[i]] <- seqsubm(channels[[i]], sm[[i]], with.missing=with.missing,
				time.varying=timeVarying, cval=cval, miss.cost=miss.cost)
		}
		## Checking correct dimension cost matrix
		else {
			if (method=="OM") {
				TraMineR.checkcost(sm[[i]], channels[[i]], with.missing=with.missing, indel=indel[i])
			} else {
				TraMineR.checkcost(sm[[i]], channels[[i]], with.missing=with.missing)
			}
			substmat_list[[i]] <- sm[[i]]
		}
		

		## Mutliply by channel weight
		substmat_list[[i]] <- cweight[i]* substmat_list[[i]]
	}
	
	## Checking that all channels have the same length
	slength1 <- seqlength(channels[[1]])
	for (i in 2:nchannels) {
		if (sum(slength1 != seqlength(channels[[i]]))>0) {
			if (!with.missing) {
				stop(" [!] some channels have sequences of different length for the same individual. Please set 'with.missing=TRUE' to nevertheless compute distances")
			} else {
				warning(" [!] some channels have sequences of different length for the same individual. Shorter sequences will be filled with missing values.")
				break
			}
		}
	}
	## ================================
	## Building the new sequence object
	## ================================
	message(" [>] building combined sequences...", appendLF=F)
	## Complex separator to ensure (hahem) unicity
	sep <- "@@@@TraMineRSep@@@@"
	maxlength=max(maxlength_list)
	newseqdata <- matrix("", nrow=numseq, ncol=maxlength)
	newseqdataNA <- matrix(TRUE, nrow=numseq, ncol=maxlength)
	for (i in 1:nchannels) {
		seqchan <- channels[[i]]
		void <- attr(seqchan, "void")
		nr <- attr(seqchan, "nr")
		for (j in 1:maxlength) {
			## No column in stslist object, filling with voids
			if (j > maxlength_list[i]) {
				newCol <- as.character(rep(void, numseq))
			}
			else {
				newCol <- as.character(seqchan[,j])
			}
			## If all channel are equal to void, then we accept as void
			newseqdataNA[,j] <- newseqdataNA[,j] & newCol == void
			## Setting void as nr
			newCol[newCol == void] <- nr
			if (i > 1) {
				newseqdata[,j] <- paste(newseqdata[,j], newCol, sep = sep)
			}
			else {
				newseqdata[,j] <- newCol
			}
		}
		
    }
	## Setting void states back to NA  (nr will be considered as a distinct state)
	newseqdata[newseqdataNA] <- NA
	
	alphabet_size <- length(unique(as.character(newseqdata))) - as.integer(sum(is.na(newseqdata))>0)
	suppressMessages(newseqdata <- seqdef(newseqdata, cpal=rep("blue", alphabet_size)))
	message(" OK")
	
	## =========================================
	## Building the new substitution cost matrix
	## =========================================
	message(" [>] computing combined substitution and indel costs...", appendLF=FALSE)
	## Build subsitution matrix and new alphabet
	alphabet <- attr(newseqdata,"alphabet")
	alphabet_size <- length(alphabet)
	## Recomputing the subsitution matrix
	if (!timeVarying) {
		newsm <- matrix(0, nrow=alphabet_size, ncol=alphabet_size)
		for (i in 1:(alphabet_size-1)) {
			statelisti <- strsplit(alphabet[i], sep)[[1]]
			for (j in (i+1):alphabet_size) {
				cost <- 0
				statelistj <- strsplit(alphabet[j], sep)[[1]]
				for (chan in 1:nchannels) {
					ipos <- match(statelisti[chan], alphabet_list[[chan]])
					jpos <- match(statelistj[chan], alphabet_list[[chan]])
					cost <- cost + substmat_list[[chan]][ipos, jpos]
				}
				newsm[i, j] <- cost
				newsm[j, i] <- cost
			}
		}
	} else {
		## Recomputing time varying substitution
		newsm <- array(0, dim=c(alphabet_size, alphabet_size, maxlength))
		for (t in 1:maxlength) {
			for (i in 1:(alphabet_size-1)) {
				statelisti <- strsplit(alphabet[i], sep)[[1]]
				for (j in (i+1):alphabet_size) {
					cost <- 0
					statelistj <- strsplit(alphabet[j], sep)[[1]]
					for (chan in 1:nchannels) {
						ipos <- match(statelisti[chan], alphabet_list[[chan]])
						jpos <- match(statelistj[chan], alphabet_list[[chan]])
						cost <- cost + substmat_list[[chan]][ipos, jpos, t]
					}
					newsm[i, j, t] <- cost
					newsm[j, i, t] <- cost
				}
			}
		}
	}
	message(" OK")
	## Indel as sum
	newindel <- sum(indel_list*cweight)
	## If we want the mean of cost..
	if (link=="mean") {
		newindel <- newindel / sum(cweight)
		newsm <- newsm / sum(cweight)
	}
	message(" [>] computing distances ...")
	## Calling seqdist
	return(seqdist(newseqdata, method=method, norm=norm, indel=newindel,
		sm=newsm, with.missing=FALSE, full.matrix=full.matrix))
	
}
