## ================================
## Convert from HSPELL to STS format
## ================================

HSPELL_to_STS <- function(seqdata, begin, end, status=NULL, fixed.status=NULL, 
	pvar=NULL, overwrite=TRUE, fillblanks=NULL, 
	tmin=NULL, tmax=NULL, id=NULL, endObs=NULL) {
	## Retrieving the list of used variable
	varlist <- c(begin, end)
	## print(varlist)
	## Process time
	if(!is.null(pvar)){
		seqdata[ , varlist] <- seqdata[ , varlist] - seqdata[ , pvar]
	}
	
	allstatus <- as.matrix(seqdata[ , status])
	## print(dim(allstatus))
	## guessing time tmin-tmax
	if(is.null(tmin)){
		tmin <- min(seqdata[ , varlist], na.rm=TRUE)
	}
	if(is.null(tmax)){
		tmax <- max(seqdata[ , varlist], na.rm=TRUE)
	}
	seqdata[, varlist] <- seqdata[ , varlist] - tmin + 1
	limit <- tmax - tmin +1
	message(" [>] sequences computed between ", tmin , " and ", tmax)
	n <- nrow(seqdata)
	seqresult <- matrix(as.character(NA), n, limit)
	seqindex<-matrix(rep(1:limit, n), nrow=n, byrow=TRUE)
	numspell <- length(begin)
	for(d in 1:numspell){
		begincolumn <- seqdata[ , begin[d]]
		endcolumn <- seqdata[ , end[d]]
		endcolumn[is.na(endcolumn)] <- tmax
		
		if(overwrite){
			cond <- !is.na(begincolumn) & seqindex>=begincolumn & seqindex<=endcolumn
			
		}else{
			cond <- !is.na(begincolumn) & seqindex>=begincolumn & seqindex<=endcolumn & is.na(seqresult)
		}
		if(!is.null(fixed.status)){
			seqresult[cond] <- fixed.status[d]
		}
		else{
			for(i in 1:nrow(seqdata)){
				seqresult[i, cond[i,]] <- allstatus[i, d]
			}
		}
	}
	if(!is.null(fillblanks)){
		seqresult[is.na(seqresult)] <- fillblanks
	}
	if(!is.null(endObs)){
		if(length(endObs)==1){
			if(is.character(endObs)){
				endObs <- seqdata[ , endObs]
			}
			else{
				endObs <- rep(endObs, n)
			}
		}
		if(!is.null(pvar)){
			endObs <- endObs - seqdata[ , pvar]
		}
		endObs <- endObs - tmin
		seqresult[seqindex >= endObs] <- NA
	}
	seqresult <- as.data.frame(seqresult)
	names(seqresult) <- paste("a", tmin:tmax, sep="")
	## setting id as rowname 
	if(!is.null(id)) {
		row.names(seqresult) <- seqdata[ , id]
	}
	
	return(seqresult)
}




