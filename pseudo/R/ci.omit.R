"ci.omit" <-
function(pseudo,tmax,causes){
	
	#calculate cum. inc. function, leave one out
	#tmax: needed only for area under CI, if specified, area is reported
	
	howmany <- nrow(pseudo)				#number of individuals
	
	ncauses <- length(causes)			#number of causes
	
	di <- matrix(NA,nrow=howmany,ncol=ncauses)
	for(jt in 1:ncauses)di[,jt] <- as.numeric(pseudo$event==causes[jt])		#indicator of  events due to certain cause

	event <- as.numeric(pseudo$event!=0)		#indicator of any event
	
	#check!!!
	td <- pseudo$time[event==1]			#times of events (any event)
	lt.temp <- c(td[-1],td[length(td)]+1)
	lt <- which(td!=lt.temp)			#in ties, use only the last value
	
	#km - i
	Y1 <- matrix(howmany:1,byrow=TRUE,ncol=howmany,nrow=howmany)
	Y2 <- matrix((howmany-1):0,byrow=TRUE,ncol=howmany,nrow=howmany)
	Y <- upper.tri(Y1,diag=FALSE)*Y1+lower.tri(Y2,diag=TRUE)*Y2
	N <- matrix(event,byrow=TRUE,ncol=howmany,nrow=howmany)
	Ndiag <- diag(diag(N))
	N <- N - Ndiag					#put zeros on diagonal
	
	cumi <- list(rep(NA,ncauses))
	for(jt in 1:ncauses){
		N1 <- matrix(di[,jt],byrow=TRUE,ncol=howmany,nrow=howmany)
		Ndiag1 <- diag(diag(N1))
		N1 <- N1 - Ndiag1			#put zeros on diagonal
		cumi[[jt]] <- N1/Y			#each line: hazard at each time, i excluded
	}
	
	kmji <- (Y-N)/Y
		
	km <- t(apply(kmji,1,cumprod))			#each line: KM estimator for theta_(-i)
	
	if(!missing(tmax)){
		tt <- matrix(pseudo$time,byrow=TRUE,nrow=nrow(pseudo),ncol=nrow(pseudo))
		#diag(tt) <- c(diag(tt[-nrow(pseudo),-1]),tmax)
		diag(tt) <- c(0,diag(tt[-1,-nrow(pseudo)]))
		tt <- tt[,pseudo$event>0,drop=FALSE]
		tt <- tt[,lt,drop=FALSE]		#pazi!!!!!!!!!!!!!!	
		tt <- cbind(tt,rep(tmax,nrow(pseudo)))
		tt <- t(apply(tt,1,diff))		#each line: times for each individual
	}
	

	#corrected value for the last time - last value carried forward 
	aje <- which(is.na(km[howmany,]))
	if(length(aje)>0){
		kir <- min(aje)
		km[howmany,kir:ncol(km)] <- km[howmany,kir-1] 
	}
	
	km <- cbind(rep(1,nrow(km)),km[,-ncol(km)])
	
	CI <- list(rep(NA,ncauses))
	
	for(jt in 1:ncauses){
	
		cit <- t(apply(cumi[[jt]]*km,1,cumsum))

		#only for deaths, one value per tie
		cit <- cit[,event==1]
		cit <- cit[,lt]
	
		if(!missing(tmax))cit <- apply(cit*tt,1,sum)
		CI[[jt]] <- cit
	}
	
	CI
}

