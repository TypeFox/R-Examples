## ====================================================
## Computing distances between sequences
## Available metrics (method):
## OM = optimal matching
## LCP = Longest Common Prefix (Elzinga)
## LCS = Longest Common Subsequence (Elzinga)
## ====================================================

seqdist <- function(seqdata, method, refseq=NULL, norm=FALSE,
	indel=1, sm=NA,	with.missing=FALSE, full.matrix=TRUE) {
	gc(FALSE)
	debut <- proc.time()

	## Checking correct arguments
	if (!inherits(seqdata,"stslist")) {
		stop(" [!] data is not a state sequence object, use 'seqdef' function to create one", call.=FALSE)
	}

    ## checking refseq
    if (inherits(refseq, "stslist")) {
        a.ref <- alphabet(refseq)
        a.seq <- alphabet(seqdata)
        if ((length(a.ref)!=length(a.seq)) | (!suppressWarnings(all(a.ref==a.seq)))) {
            stop(" [!] same alphabet must be assigned to 'refseq' and 'seqdata'!", call.=FALSE)
        }
    }

	if (method=="OMopt") {
		method <- "OM"
		optimized <- TRUE
	}
	else if (method=="LCSopt") {
		method <- "LCS"
		optimized <- TRUE
	}
	else {
		optimized <- FALSE
	}
	metlist <- c("OM","LCP", "LCS", "RLCP", "DHD", "HAM")
	if (missing(method)) {
		stop(" [!] You should specify a method to compute the distances. It must be one of: ", paste(metlist,collapse=" "))
	}
	if (!method %in% metlist) {
		stop(" [!] Method must be one of: ", paste(metlist,collapse=" "), call.=FALSE)
	}
	## Taking care of correct normalization settings
	if (is.logical(norm)) {
		if (norm) {
		## Normalize using Elzinga for LCP, LCS, RLCP
			if (method %in% c("LCP", "LCS", "RLCP")) {
				norm <- 2
			} else {## Normalize using Abbott for OM, HAM, and DHD
				norm <- 1
			}
      	} else {
      		norm <- 0
		}
	} else if (is.character(norm)) {
		## Using normalization name
		## Cast to integer for c code and normdist function
		## Match return the position, removing 1 to start at zero
		normIndice <- match(norm, c("none", "maxlength", "gmean", "maxdist", "YujianBo")) -1

		if (is.na(normIndice)) {  ##Not found
			stop(" [!] Unknown distance normalization method ", norm)
		}
		norm <- normIndice
	} else {
		stop(" [!] Unknown distance normalization method ", norm)
	}
	## Checking missing values
	if (!with.missing && any(seqdata==attr(seqdata,"nr"))) {
		stop("found missing values in sequences, please set 'with.missing=TRUE' to nevertheless compute distances")
	}
	if (method == "OM" && is.character(sm)) {
		if (sm == "TRATE" || sm =="CONSTANT") {
			sm <- seqsubm(seqdata, method=sm, cval=2, with.missing=with.missing, miss.cost=2)
		} else {
			stop(" [!] Unknown method ", sm, " to compute substitution costs")
		}
	}
	if (method == "OM" && is.na(sm)) {
		stop(" [!] sm=NA not allowed with method='OM'. Specify a valid sm value!")
    }

	
	## =====================
	## Base information
	## =====================
	
	n <- nrow(seqdata)
	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)
	message(" [>] ",n," sequences with ", alphsize, " distinct events/states")
	## Gaps in sequences
	if (with.missing) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
		alphsize <- length(alphabet)
		message(" [>] including missing value as additional state" )
	}
	
	## Checking methods that are treated the same
	methodname <- method
	if (method == "LCS") {
		method <- "OM"
		sm <- suppressMessages(seqsubm(seqdata, method="CONSTANT", cval=2, with.missing=with.missing, miss.cost=2))
		indel <- 1
	} else if (method == "HAM") {
		method <- "DHD"
		
		if (!is.null(dim(sm))) {
			TraMineR.checkcost(sma=sm, seqdata=seqdata, with.missing=with.missing)
			if(is.matrix(sm)){
				costs <- array(0, dim=c(alphsize, alphsize, ncol(seqdata)))
				for(i in 1:ncol(seqdata)){
					costs[,,i] <- sm
				}
				sm <- costs
			}
		} else {
			sm <- suppressMessages(seqsubm(seqdata, "CONSTANT", cval=1, with.missing=with.missing,
				miss.cost=1, time.varying=TRUE))
		}
	}


	## ===========================
	## Checking correct size of sm
	## ===========================
	
	## Checking if substitution cost matrix contains values for each state
	## and if the triangle inequality is respected
	if (method=="OM") {
		TraMineR.checkcost(sma=sm, seqdata=seqdata, with.missing=with.missing, indel=indel)
	}
	## Checking if substitution cost matrix contains values for each state
	## and if the triangle inequality is respected
	if(methodname == "DHD"){
		## User entered substitution cost
		if(length(sm)>1){
			TraMineR.checkcost(sma=sm, seqdata=seqdata, with.missing=with.missing)
		}
		else {
			sm <- seqsubm(seqdata, "TRATE", cval=4, with.missing=with.missing,
					miss.cost=4, time.varying=TRUE)
		}
	}


	## ==============
	## Preparing data
	## ==============
	seqdata <- seqnum(seqdata, with.missing=with.missing)

	## Selecting distinct sequences only and saving the indexes
	dseq <- unique(seqdata)
	mcorr <- match(seqconc(seqdata), seqconc(dseq))

	nd <- nrow(dseq)
	message(" [>] ", nd," distinct sequences")

	slength <- seqlength(dseq)

	dseq <- seqasnum(dseq, with.missing=with.missing)

	message(" [>] min/max sequence length: ",min(slength),"/",max(slength))

	## ===================================
	## Preparing data for Hamming distance
	## ===================================
	if (method=="DHD") {
		if(length(unique(slength))>1) {
			## Hamming is not defined for sequence of different length
			stop(methodname, " distance can only be computed between sequences of equal length")
		}
		## Here we use the indel parameter for the C function
		## But it contain the maximum possible cost of the hamming distance
		indel <- 0
		for (i in 1:max(slength)) {
			indel <- indel + max(sm[,,i])
		}
	}
	message(" [>] computing distances using ", methodname,
		ifelse(norm!=0," normalized", ""), " metric")
	
	## Function and arguments
	if (!missing(refseq) && !is.null(refseq)) {
		distances <- TraMineR.seqdist.refseq(seqdata, method, refseq,
			norm, indel, sm, alphsize, nd, dseq, slength, mcorr, with.missing)	
	}
	else {
		distances <- TraMineR.seqdist.all(seqdata, method,
			norm, indel, sm, alphsize, nd, dseq, slength, mcorr, optimized)
	}

	fin <- proc.time()
	totaltime <- format(round(difftime(as.POSIXct(sum(fin[1:2]), origin="1960-01-01"), as.POSIXct(sum(debut[1:2]), origin="1960-01-01")), 3))
	message(" [>] total time: ", totaltime)

	if (full.matrix && inherits(distances, "dist")) {
		return(dist2matrix(distances))
	}
	else {
		return(distances)
	}
}
