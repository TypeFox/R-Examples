"pseudomean" <-
function(time, event, tmax){

	if(any(is.na(time)))
		stop("missing values in 'time' vector")
		
	if(any(time<0))
		stop("'time' must be nonnegative")
	
	if(any(is.na(event)))
		stop("missing values in 'event' vector")
	
	if(any(event!=0 & event!=1))
		stop("'event' must be a 0/1 variable (alive/dead)")
		
	if(missing(tmax)) 
		tmax <- max(time[event==1])

	if(length(tmax)>1)
		stop("Only one value should be specified for 'tmax' - the maximum time")

	if(is.na(tmax))
		stop("missing value of 'tmax'")

	
	rmtime <- ifelse(time >= tmax ,tmax,time)			#if times are greater than tmax, set them to tmax
	rmdead <- ifelse(time > tmax ,0, event)			#if times are greater than tmax, set censoring status to 0

	if(sum(rmdead)==0)
		stop("no events occured before time 'tmax'")
    
	howmany <- length(rmtime)
    	
    	## preparing the data
  	pseudo <- data.frame(id=1:howmany,time=rmtime,event=rmdead)
  
  	# sort in time, if tied, put events before censoring
  	pseudo <- pseudo[order(pseudo$time,-pseudo$event),]

    	# RM, leave one out
	RM.omit <- surv.omit(pseudo,tmax)
	
	#RM, all cases
	RM.tot <- surv.tot(pseudo,tmax)

	# pseudo-observations
	pseu <- howmany*RM.tot - (howmany-1)*RM.omit

	#back to original order
	pseu <- pseu[order(pseudo$id)]

	return(pseu)    
}

