seqstart <- function(seqdata, data.start, new.start, tmin=NULL, tmax=NULL, missing=NA){
	warning(" [!] Please check your results. This function needs further testing.")
	new.index <- as.integer(data.start - new.start+1)
	if(length(new.index)!=nrow(seqdata)){
		if(length(new.index)==1){
			new.index <- rep(new.index, nrow(seqdata))
		}
		else{
			stop(" [!] number of individual in seqdata, data.start and/or new.start mismatch.")
		}
	}
	if(any(new.index <1)){
		correction <- -min(new.index)+1
		new.index <- new.index + correction
	} else {
		correction <- 0
	}
	if(is.null(tmin)){
		tmin <- max(min(new.index), 1)
	}
	tmin <- tmin+correction
	if(is.null(tmax)){
		tmax <- max(new.index)+ncol(seqdata) -1
	}
	tmax <- tmax +correction
	if(tmax<0){
		stop("[!] There are no data in the specified new time frame.")
	}
	## cat(tmin, tmax, correction)
	##new.index.mat <- new.index - tmin + 1
	new.data <- matrix(as.character(missing), ncol=(tmax-tmin+1), nrow=nrow(seqdata))
	seqdatadim <- dim(seqdata)
	seqdata <- as.character(as.matrix(seqdata))
	dim(seqdata) <- seqdatadim
	#tmin <- tmin + 1 - firstyear
	## print(tmin:tmax-correction-1)
	return(.Call(TMREXTRAS_tmrextrasseqstart, seqdata, new.data, as.integer(new.index-tmin)))
}
