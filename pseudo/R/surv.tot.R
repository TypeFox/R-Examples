"surv.tot" <-
function(pseudo,tmax){

	# calculate Kaplan - Meier, all cases

	howmany <- nrow(pseudo)
	
	td <- pseudo$time[pseudo$event==1]
	lt.temp <- c(td[-1],td[length(td)]+1)
	lt <- which(td!=lt.temp)
	
	#km - i
	Y <- howmany:1
	N <- pseudo$event
	
	kmji <- (Y-N)/Y
		
	km <- cumprod(kmji)
	
	if(!missing(tmax)){
		tt <- pseudo$time[pseudo$event==1]
		tt <- tt[lt]
		tt <- c(0,tt,tmax)
		tt <- diff(tt)
	}
	
	#only for deaths, one value per tie
	km <- km[pseudo$event==1]
	km <- km[lt]
	if(!missing(tmax)){
		km <- sum(c(1,km)*tt)
	}
	km	
}

