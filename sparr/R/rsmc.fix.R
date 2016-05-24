
rsmc.fix <- function(object,hpsim,ITER,corrgrid,comment){
	datapoolx <- append(as.vector(object$f$data[,1]),as.vector(object$g$data[,1]))
	datapoolx <- rep(datapoolx,append(object$f$counts,object$g$counts))
	datapooly <- append(as.vector(object$f$data[,2]),as.vector(object$g$data[,2]))
	datapooly <- rep(datapooly,append(object$f$counts,object$g$counts))
	
	n1 <- sum(object$f$counts)
	n2 <- sum(object$g$counts)
	datapool <- data.frame(cbind(datapoolx,datapooly))
	range.n <- 1:nrow(datapool)
	
	res <- sqrt(length(corrgrid))
	WIN <- object$f$WIN
	rmat <- matrix(as.vector(t(object$rsM))[corrgrid],res,res,byrow=T)
	mcmat_upper <- matrix(1,res,res)
		
	if(comment) pb <- txtProgressBar(0,ITER-1,style=3)
	for(i in 1:(ITER-1)){
		casetemp <- sample(range.n,n1)
		contemp <- range.n[-casetemp]
		
		if(is.null(hpsim)){
			hpf <- object$f$pilotH
			hpg <- object$g$pilotH
			
		} else if(is.function(hpsim)){
			hp <- hpsim(datapool[casetemp,],datapool[contemp,])

			if(!is.numeric(hp)) stop("if a function, 'hpsim' must return a numeric vector of length 1 or 2")
			
			if(length(hp)==1){
				hpf <- hpg <- hp
			} else if(length(hp)==2){
				hpf <- hp[1]
				hpg <- hp[2]
			} else {
				stop("if a function, 'hpsim' must return a numeric vector of length 1 or 2")
			}
						
		} else if(is.numeric(hpsim)){
			if(length(hpsim)!=2 && length(hpsim)!=1) stop("if numeric, 'hpsim' must be a vector of length 1 or 2")

			if(length(hpsim)==1){
				hpf <- hpg <- hpsim
			} else {
				hpf <- hpsim[1]
				hpg <- hpsim[2]
			}

		} else {
			stop("invalid 'hpsim' argument")
		}		
		
		casedens <- bivariate.density(data=datapool[casetemp,],pilotH=hpf,adaptive=F,res=res,WIN=WIN,comment=F)
		condens <- bivariate.density(data=datapool[contemp,],pilotH=hpg,adaptive=F,res=res,WIN=WIN,comment=F)
		
		risktemp <- risk(casedens,condens,log=object$log,plotit=F)$rsM
			
		mcmat_upper <- mcmat_upper+(risktemp>=rmat)
		if(comment) setTxtProgressBar(pb,i)
	}
	if(comment) close(pb)

	return(mcmat_upper/ITER)
}
