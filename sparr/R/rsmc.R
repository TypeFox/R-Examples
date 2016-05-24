rsmc <- function(object,pooled,ITER,corrgrid,comment){
	datapoolx <- append(as.vector(object$f$data[,1]),as.vector(object$g$data[,1]))
	datapoolx <- rep(datapoolx,append(object$f$counts,object$g$counts))
	datapooly <- append(as.vector(object$f$data[,2]),as.vector(object$g$data[,2]))
	datapooly <- rep(datapooly,append(object$f$counts,object$g$counts))
	
	n1 <- sum(object$f$counts)
	n2 <- sum(object$g$counts)
	datapool <- data.frame(cbind(datapoolx,datapooly))
	range.n <- 1:nrow(datapool)
	
	pilotH_f <- object$f$pilotH
	pilotH_g <- object$g$pilotH
	if(is.na(object$f$globalH)) globalH_f <- pilotH_f
	else globalH_f <- object$f$globalH
	if(is.na(object$g$globalH)) globalH_g <- pilotH_g
	else globalH_g <- object$g$globalH
	adaptive <- (length(object$f$hypoH)>1)
	res <- sqrt(length(corrgrid))
	WIN <- object$f$WIN
	rmat <- matrix(as.vector(t(object$rsM))[corrgrid],res,res,byrow=T)
	mcmat_upper <- matrix(1,res,res)
	#mcmat_lower <- matrix(1,res,res)
	
	
	if(comment) cat("\nMonte-Carlo iteration no.\n")
	for(i in 1:(ITER-1)){
		if(comment){
			if(i%%10==0) cat(i," ")
			if(i%%100==0) cat("\n")
		}
		
		casetemp <- sample(range.n,n1)
		contemp <- range.n[-casetemp]
		
		casedens <- bivariate.density(data=datapool[casetemp,],pilotH=pilotH_f,globalH=globalH_f,adaptive=adaptive,res=res,WIN=WIN,gamma=pooled$gamma,comment=F)
		condens <- bivariate.density(data=datapool[contemp,],pilotH=pilotH_g,globalH=globalH_g,adaptive=adaptive,res=res,WIN=WIN,gamma=pooled$gamma,comment=F)
		
		risktemp <-risk(casedens,condens,log=object$log,plotit=F)$rsM
			
		mcmat_upper <- mcmat_upper+(risktemp>=rmat)
		#mcmat_lower <- mcmat_lower+(risktemp<=rmat)
	}

	return(mcmat_upper/ITER)
}