## ================================
## Convert from STS to TSE format
## Version 1.3
## ================================

STS_to_TSE <- function(seqdata, id=NULL, tevent) {
	nseq <- nrow(seqdata)
	slength <- ncol(seqdata)

	## Ceci accélère énormément la fonction
	seqdata <- as.matrix(seqdata)

	message(" [>] converting ",nseq, " sequences to TSE format ...")
 
	if (is.null(id)) id <- 1:nseq
	maxsize <- nseq*slength
	ids <- vector(mode=mode(id),length=maxsize)
	times <- numeric(length=maxsize)
	events <- vector(mode=mode(tevent),length=maxsize)
	myi <- 1

	levent <- list()
	nl <- nrow(tevent)
	nc <- ncol(tevent)

	## Storing the cells of tevent in a list
	for (i in 1:nc){
		for (j in 1:nl) {
			if (is.character(tevent[j,i])){
				ll <- strsplit(tevent[j,i],",")[[1]]
				if(length(ll)==0){
					ll <- NA
				}
				levent[[paste(rownames(tevent)[j],">",colnames(tevent)[i],sep="")]] <- ll
			} else {
				levent[[paste(rownames(tevent)[j],">",colnames(tevent)[i],sep="")]] <- as.character(tevent[j,i])
			}
		}
	}
	for (i in 1:nseq) {
		## First status=> entrance event (diagonal of tevent)
        s1 <- seqdata[i,1]
		e1 <- levent[[paste(s1,">",s1,sep="")]]

		if (!is.null(e1) && !is.na(e1[1])) { ## if NA, we don't generate an event
			for (k in e1){
				ids[myi] <- id[i]
				times[myi] <- 0
				events[myi] <- k
				myi <- myi+1
			}
	     } #end if
		## Rest of the sequence
		for (j in 1:(slength-1)) {
			s1 <- seqdata[i,j] ## Status at time t
			s2 <- seqdata[i,j+1] ## Status at time t+1

			if (!is.na(s1) && !is.na(s2) && s2!=s1 && s1!="" && s2!="") {
				ev <- levent[[paste(s1,">",s2,sep="")]]
				if(!is.na(ev[1])){
					for (k in ev){
						ids[myi] <- id[i]
						times[myi] <- j
						events[myi] <- k
						myi <- myi+1
					}
				}
			}
		}
	}

	sel <- 1:(myi-1)
	trans <- data.frame(id=ids[sel],time=times[sel],event=events[sel])
	return(trans)
}


