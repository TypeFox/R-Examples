"surv.omit" <-
function(pseudo,tmax){
	
	# calculate Kaplan - Meier, leave one out 
	howmany <- nrow(pseudo)
	
	td <- pseudo$time[pseudo$event==1]
	lt.temp <- c(td[-1],td[length(td)]+1)
	lt <- which(td!=lt.temp)
	
	#km - i
	Y1 <- matrix(howmany:1,byrow=TRUE,ncol=howmany,nrow=howmany)
	Y2 <- matrix((howmany-1):0,byrow=TRUE,ncol=howmany,nrow=howmany)
	Y <- upper.tri(Y1,diag=FALSE)*Y1+lower.tri(Y2,diag=TRUE)*Y2
	N <- matrix(pseudo$event,byrow=TRUE,ncol=howmany,nrow=howmany)
	Ndiag <- diag(diag(N))
	N <- N - Ndiag
	
	kmji <- (Y-N)/Y
		
	km <- t(apply(kmji,1,cumprod))
	
	if(!missing(tmax)){
		tt <- matrix(pseudo$time,byrow=TRUE,nrow=nrow(pseudo),ncol=nrow(pseudo))
		#diag(tt) <- c(diag(tt[-nrow(pseudo),-1]),tmax)
		diag(tt) <- c(0,diag(tt[-1,-nrow(pseudo)]))
		tt <- tt[,pseudo$event==1,drop=FALSE]
		tt <- tt[,lt,drop=FALSE]
		tt <- cbind(rep(0,nrow(pseudo)),tt,rep(tmax,nrow(pseudo)))
		tt <- t(apply(tt,1,diff))
	}
	
	
	#corrected value for the last time - last value carried forward 
	aje <- which(is.na(km[howmany,]))
	if(length(aje)>0){
		kir <- min(aje)
		km[howmany,kir:ncol(km)] <- km[howmany,kir-1] 
	}
	
	#only for deaths, one value per tie
	km <- km[,pseudo$event==1,drop=FALSE]
	km <- km[,lt,drop=FALSE]
	if(!missing(tmax)){
		km <- apply(cbind(rep(1,nrow(pseudo)),km)*tt,1,sum)
	}
	km	
}

