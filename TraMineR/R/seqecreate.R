## ========================================
## Create events objects
## ========================================


#tmrsequence<-function(id,timestamp,event){
#  .Call("tmrsequence",as.integer(id),as.double(timestamp),as.integer(event), PACKAGE="TraMineR")
#}
seqecreate <- function(data=NULL, id=NULL,timestamp=NULL, event=NULL, endEvent=NULL,
						tevent="transition", use.labels=TRUE, weighted=TRUE){
	return(seqecreate.internal(data=data, id=id, timestamp=timestamp, event=event,
								endEvent=endEvent, tevent=tevent, use.labels=use.labels,
								order.before=FALSE, weighted=weighted))
}
seqecreate.internal <- function(data, id, timestamp, event, endEvent, tevent,
								use.labels, order.before, weighted){
	if (!is.null(data)) {
		if (inherits(data,"stslist")) {
			if (!is.matrix(tevent)) {
				if (is.character(tevent)) {
					tevent <- seqetm(data, method=tevent, use.labels=use.labels)
				}else{
					tevent <- seqetm(data, use.labels=use.labels)
				}
			}
			data.tse <- suppressMessages(seqformat(data, from='STS',to='TSE', tevent=tevent))
			id <- data.tse$id
			timestamp <- data.tse$time
			event <- data.tse$event
		} else if(is.data.frame(data)) {
			dname <- names(data)
			if("id" %in% dname && ("timestamp" %in% dname ||"time" %in% dname)&& "event" %in% dname){
				id <- data[ ,"id"]
				event <- data[ ,"event"]
				if ("timestamp" %in% dname) {
					timestamp <- data[ ,"timestamp"]
				} else {
					timestamp <- data[ ,"time"]
				}
			}
		}
	}
	if (is.null(id)) {
		stop(" [!] Could not find an id argument")
	}
	if(is.null(timestamp))	{
		stop(" [!] Could not find a timestamp argument")
	}
	if(is.null(event)) {
		stop(" [!] Could not find an event argument")
	}
	if (any(is.na(id)) || any(is.na(timestamp)) || any(is.na(event))) {
		stop(" [!] Missing values not supported")
	}
	#  warning("Event sequence analysis module is still experimental", call.=FALSE)
	classname <- c("seqe")
	intEvent <- NULL
	if(!is.factor(event)){
		event <- factor(event)
	}
	if (is.factor(event)) {
		dictionnary <- levels(event)
		if (!is.null(endEvent)) {
			for(i in 1:length(dictionnary)){
				if (dictionnary[i] == endEvent) {
					intEvent <- i
				}
			}
			if (is.null(intEvent)) {
				stop(" [!] endEvent not found in event dictionary")
				return(invisible())
			}
		}
	} else {
		dictionnary <- c()
	}
	if(any(grepl("\\(|\\)|,", dictionnary))){
		warning(" [!] some of your events contain '(', ')' or ',' characters. The search of specific subsequences may not work properly.")
	}
	id <- as.integer(id)
	timestamp <- as.double(timestamp)
	event <- as.integer(event)

	
	ret <- .Call(TMR_tmrsequenceseveral, as.integer(id),
		as.double(timestamp), as.integer(event),
		as.integer(c(intEvent)), classname, as.character(dictionnary))

	class(ret) <- c("seqelist", "list")
	if(inherits(data, "stslist")){
		seqelength(ret) <- seqlength(data)
		ww <- attr(data, "weights")
		if(!is.null(ww) && weighted){
			seqeweight(ret) <- ww
		}
	}
	if(length(ret) != length(unique(id))){
		stop(" [!] Events not grouped by id! See seqecreate help page.")
	}
	return(ret)
}
#SEXP tmrsequence(SEXP idpers, SEXP time, SEXP event, SEXP classname, SEXP seq)

seqecreatesub <- function(subseq, seqe){
#  warning("Event sequence analysis module is still experimental", call.=FALSE)
	if (!is.seqelist(seqe)) {
		stop(" [!] seqe should be a seqelist. See help on seqecreate.")
	}
	classname <- c("seqe")
	irow <- 1
	ret <- list()
	codebase <- levels(seqe)
	iseq <- 1
	for (subseqstr in subseq) {
		mystr <- gsub("(^\\()|(\\)$)", "", unlist(strsplit(subseqstr, "\\)[[:space:]]*-[[:space:]]*\\(")))
		timestamp <- numeric()
		events <- integer()
		tindex <- 1
		irow <- 1
		for (m in mystr) {
			mm <- unlist(strsplit(m, "\\,"))
			for (mmm in mm) {
				mmm <- sub('[[:space:]]+$', '', mmm)
				mmm <- sub('^[[:space:]]+', '', mmm)
				ecode <- charmatch(mmm, codebase)
				if (is.na(ecode)) {
					stop(" [!] Couldn't interpret '", mmm,"' as an event. It should be in (", paste(codebase, collapse=","),")")
				}
				timestamp[irow] <- tindex
				events[irow] <- ecode
				irow <- 1 + irow
			}
			tindex <- tindex+1
		}
		timestamp <- as.double(timestamp)
		events <- as.integer(events)
		sortedindex <- order(timestamp, events)
		
		ret[[iseq]]<-.Call(TMR_tmrsequence, as.integer(-1),
			as.double(timestamp[sortedindex]), as.integer(events[sortedindex]),
			classname, seqe[[1]])
			iseq <- iseq + 1
	}
#  e<-factor(event,levels=levels(seq))
 # ret<-list(.Call("tmrsequence",as.integer(-1),as.double(timestamp),as.integer(e),classname,seq, PACKAGE="TraMineR"))
	class(ret) <- c("seqelist", "list")
	return(ret)
}
