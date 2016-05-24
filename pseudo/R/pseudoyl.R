"pseudoyl" <-
function(time, event, tmax){

	if(any(is.na(time)))
		stop("missing values in 'time' vector")
		
	if(any(time<0))
		stop("'time' must be nonnegative")
	
	if(any(is.na(event)))
		stop("missing values in 'event' vector")
	
	if(!is.integer(event[event!=0])| any(event<0))
		stop("Events should be denoted by integer values, 0 is used for censoring")

	if(missing(tmax)) 
		tmax <- max(time[event!=0])
	
	if(length(tmax)>1)
		stop("Only one value should be specified for 'tmax' - the maximum time")
	
	if(is.na(tmax))
		stop("missing value of 'tmax'")

	
	rmtime <- ifelse(time >= tmax ,tmax,time)			#if times are greater than tmax, set them to tmax
	rmevent <- ifelse(time > tmax ,0, event)			#if times are greater than tmax, set censoring status to 0

	if(sum(rmevent)==0)
		stop("no events occured before time 'tmax'")
    
	howmany <- length(rmtime)

	causes <- sort(unique(event[event!=0]))	#distinct codes for causes (censorings excluded)
	ncauses <- length(causes)			#number of causes
	
	pseudo <- data.frame(id=1:howmany,time=rmtime,event=rmevent) 	#form a data frame, use new times and censoring statuses
	# sort in time, put events before censorings (if happening simultaneously)
	pseudo <- pseudo[order(pseudo$time, -pseudo$event), ]
    		 
    
	# years lost, leave one out
	L.omit <- ci.omit(pseudo,tmax,causes=causes)
	
	# years lost, all cases
	L.tot <- ci.tot(pseudo,tmax,causes=causes)

	# pseudo-observations
	out <- NULL
	out$cause <- causes
	out$pseudo <- list(rep(NA,ncauses))
	for(jt in 1:ncauses){
		ps <- howmany*L.tot[[jt]] - (howmany-1)*L.omit[[jt]]
		out$pseudo[[jt]] <- ps[order(pseudo$id)] 		#back to original order
	}
	names(out$pseudo) <- paste("cause",causes,sep="")
	return(out)    
}

