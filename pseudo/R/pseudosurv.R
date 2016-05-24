"pseudosurv" <-
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
		tmax <- unique(time[event==1])
	
	if(any(is.na(tmax)))
		stop("missing values in 'tmax' vector")
	
	if (sum(tmax > max(time)) > 0) 
	   stop ("tmax greater than largest observation time")
	   
	tmax <- sort(tmax)
	ltmax <- length(tmax)
	howmany <- length(time)
	

	## preparing the data
	pseudo <- data.frame(id=1:howmany,time=time,event=event)

	# sort in time, if tied, put events before censoring
	pseudo <- pseudo[order(pseudo$time,-pseudo$event),]

	# time points chosen	
	tu <- unique(pseudo$time[pseudo$event==1])
	ltu <- length(tu)
	tu <- matrix(tu,byrow=TRUE,ncol=ltu,nrow=ltmax)
	tgiven <- matrix(tmax,byrow=FALSE,ncol=ltu,nrow=ltmax)
	inx <- apply(tgiven>=tu,1,sum)
	
	# KM, leave one out
	KM.omit <- surv.omit(pseudo)
	KM.omit <- KM.omit[,inx,drop=FALSE]
	
	# KM, all cases
	KM.tot <- surv.tot(pseudo)[inx]
	KM.tot <- matrix(KM.tot,byrow=TRUE,nrow=howmany,ncol=length(tmax))
	
	# pseudo-observations
	pseu <- howmany*KM.tot - (howmany-1)*KM.omit
	
	pseu <- as.data.frame(pseu)
	row.names(pseu) <- pseudo$id
	names(pseu) <- 	paste("time",tmax,sep=".")
	
	# back to original order
	out <- NULL
	out$time <- tmax
	out$pseudo <- as.matrix(pseu[order(pseudo$id),])		#back to original order
	return(out)

}

