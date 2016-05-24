"pseudoci" <-
function(time, event, tmax){

	if(any(is.na(time)))
		stop("missing values in 'time' vector")
		
	if(any(time<0))
		stop("'time' must be nonnegative")
	
	if(any(is.na(event)))
		stop("missing values in 'event' vector")
	
	if(any(event<0))
		stop("Events should be denoted by integer values, 0 is used for censoring")
		
	if(missing(tmax)) 
		tmax <- unique(time[event!=0])
	
	if(any(is.na(tmax)))
		stop("missing values in 'tmax' vector")
	
	if (sum(tmax > max(time)) > 0) 
	   stop ("tmax greater than largest observation time")
   
	tmax <- sort(tmax)				#times for evaluation
	ltmax <- length(tmax)				#number of times for evaluation
	howmany <- length(time)				#number of individuals
	
	causes <- sort(unique(event[event!=0]))	#distinct codes for causes (censorings excluded)
	ncauses <- length(causes)			#number of causes
	
	pseudo <- data.frame(id=1:howmany,time=time,event=event)
	# sort in time, put events before censorings (if happening simultaneously)
	pseudo <- pseudo[order(pseudo$time, -pseudo$event), ]

	
	# time points chosen	
	tu <- unique(pseudo$time[pseudo$event!=0])			#unique time of failures
	ltu <- length(tu)						#number of unique times of failures
	tu <- matrix(tu,byrow=TRUE,ncol=ltu,nrow=ltmax)			#matrix: unique times of failures in each row, (each evaluation time) x (each failure time) 
	tgiven <- matrix(tmax,byrow=FALSE,ncol=ltu,nrow=ltmax)		#matrix: evaluation times in columns, (each evaluation time) x (each failure time)  
	inx <- apply(tgiven>=tu,1,sum)					#for each evaluation time: the number of events before it

	# CI, leave one out
	pseu <- ci.omit(pseudo,causes=causes)						#calculate the values of CI functions for 'leave one out'
	pseut <- ci.tot(pseudo,causes=causes)						#calculate the values of CI for the whole sample
	
	out <- NULL
	out$time <- tmax
	out$cause <- causes
	out$pseudo <- list(rep(NA,ncauses))
	
	for(jt in 1:ncauses){						#for each cause
		ci <- pseu[[jt]][,inx,drop=FALSE]				#only at times for evaluation								
		citot <- matrix(pseut[[jt]][inx],byrow=TRUE,ncol=ncol(ci),nrow=nrow(ci))	#only at time for evaluation, put in the same format as the 'leave one out'
		ps <- as.data.frame(howmany*citot - (howmany-1)*ci)			#pseudo values calculation	
		row.names(ps) <- pseudo$id
		names(ps) <- paste("time",tmax,sep=".")
		out$pseudo[[jt]] <- as.matrix(ps[order(pseudo$id),])		#back to original order
		
		
	}
	names(out$pseudo) <- paste("cause",causes,sep="")
	out
}

