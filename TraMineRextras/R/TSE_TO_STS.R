TSE_to_STS <- function(seqdata, id=1, timestamp=2, event=3, stm=NULL, tmin=1, tmax=NULL, firstState="None") {
	
	event <- as.character(seqdata[, event])
	id <- as.character(seqdata[, id])
	timestamp <- as.numeric(seqdata[, timestamp])+1
	nseq <- length(unique(id))
	uid <- character(nseq)
	if(is.null(stm)){
		stm <- seqe2stm(unique(event))
	}
	
	#eorder <- order(id, timestamp, event)
	#event <- event[eorder]
	#timestamp <- timestamp[eorder]+1
	#id <- id[eorder]
	#uid <- unique(id)	
	if(is.null(tmax)){
		stop(" [!] You should use the tmax arguments to specify the length of the sequence") 
	}
	sts <- matrix(firstState, nrow=nseq, ncol=tmax)
	oldID <- NULL
	iID <- 1
	myi <- 1
	nid <- length(id)
	while(iID <= nid) {
		oldID <- id[iID]
		while(iID <= nid && id[iID]==oldID){
			iID <- iID + 1
		}
		if(iID==nid){
			break;
		}
		mid <- oldID
		uid[myi] <- mid
		cond <- mid==id
		mevent <- event[cond]
		mtime <- timestamp[cond]
		oo <- order(mtime, mevent)
		mevent <- mevent[oo]
		mtime <- mtime[oo]
		prevState <- firstState
		prevTime <- 1
		for(i in 1:length(mevent)){
			if (mevent[i] %in% colnames(stm)) {
		  #		print(mevent[i])
				tt <- min(tmax, floor(mtime[i]))
				sts[myi, prevTime:tt] <- prevState
		#		message("Indice, ", prevTime:tt, " prevState ", prevState)
				prevTime <- tt
				prevState <- stm[prevState, mevent[i]]
			}
		}
		sts[myi, prevTime:tmax] <- prevState
		#message("Indice, ", prevTime:tmax, " prevState ", prevState)
		#print(sts[myi,])
		myi <- myi +1
	}
	
	sts <- sts[,tmin:tmax]
	rownames(sts) <- uid
	colnames(sts) <- paste("a", tmin:tmax, sep="")
	return(as.data.frame(sts))
}
