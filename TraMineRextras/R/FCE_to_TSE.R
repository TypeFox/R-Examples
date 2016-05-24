FCE_to_TSE <- function(seqdata, id=NULL, cols, eventlist=NULL, firstEvent=NULL){
	addbirth <- !is.null(firstEvent)
	if (is.null(eventlist)) {
		if (is.numeric(cols)) {
			eventlist <- names(seqdata)[cols]
		}
		else {
			eventlist <- cols
		}
	}
	num_rows <- sum(!is.na(seqdata[,cols]))
	if (addbirth) {
		num_rows <- num_rows + nrow(seqdata)
	}
	if (is.null(id)) id <- 1:nrow(seqdata)
	if(length(id)==1) id <- seqdata[, id]
	ids <- vector(mode=mode(id),length=num_rows)
	times <- numeric(length=num_rows)
	events <- character(length=num_rows)
	myi <- 1
	for (i in 1:nrow(seqdata)) {
		firsti <- myi
		if (addbirth) {
			times[myi] <- 0
			events[myi] <- firstEvent
			myi <- myi+1
		}
		for(k in 1:length(cols)) {
			if(!is.na(seqdata[i,cols[k]])) {
				times[myi] <- seqdata[i,cols[k]]
				events[myi] <- eventlist[k]
				myi <- myi+1
			}
		}
		## reorder for ascending time/event
		sel <- firsti:(myi-1)
		ids[sel] <- id[i]
		timessel <- times[sel]
		eventsel <- events[sel]
		oo <- order(timessel, eventsel)
		times[sel] <- timessel[oo]
		events[sel] <- eventsel[oo]
	}
	sel <- 1:(myi-1)
	trans <- data.frame(id=ids[sel],time=times[sel],event=events[sel])
	return(trans)
}