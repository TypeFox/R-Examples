"ci.tot" <-
function(pseudo,tmax,causes){
	#calculate cum. inc. function, all cases
	#current ties handling: if two individuals fail at the same time, we split ties (use original ordering)
	#tmax: needed only for area under CI, if specified, area is reported
	
	
	howmany <- nrow(pseudo)				#number of individuals
	
	ncauses <- length(causes)			#number of causes
		
	di <- matrix(NA,nrow=howmany,ncol=ncauses)
	for(jt in 1:ncauses)di[,jt] <- as.numeric(pseudo$event==causes[jt])		#indicator of  events due to certain cause

	event <- as.numeric(pseudo$event!=0)		#indicator of any event	
	
	td <- pseudo$time[event==1]			#times of events (any event)
	lt.temp <- c(td[-1],td[length(td)]+1)
	lt <- which(td!=lt.temp)			#in ties, use only the last value

	#km - i
	Y <- howmany:1					
	N <- event
	
	cumi <- list(rep(NA,ncauses))
	for(jt in 1:ncauses)cumi[[jt]] <- di[,jt]/Y			#each line: hazard at each time, i excluded
	
	kmji <- (Y-N)/Y					#1- dN/Y
		
	km <- cumprod(kmji)				#Kaplan-Meier at each time 
	
	km <- c(1,km[-length(km)])			#KM at t- (just before that time)
	
	if(!missing(tmax)){				
		tt <- pseudo$time[pseudo$event>0]
		tt <- tt[lt]
		tt <- c(tt,tmax)
		tt <- diff(tt)				#difference between times
	}

	CI <- list(rep(NA,ncauses))
	
	for(jt in 1:ncauses){
	
		cit <- cumsum(cumi[[jt]]*km)				#cuminc: cumsum(S(t-)*hazard)

		#only for deaths, one value per tie
		cit <- cit[event==1]
		cit <- cit[lt]
	
		if(!missing(tmax))cit <- sum(cit*tt)
		CI[[jt]] <- cit
	}
CI
}

