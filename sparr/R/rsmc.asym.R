
rsmc.asym <- function(object,hpsim,h0sim,ITER,corrgrid,comment){ 
	datapoolx <- c(as.vector(object$f$data[,1]),as.vector(object$g$data[,1]))
	datapoolx <- rep(datapoolx,c(object$f$counts,object$g$counts))
	datapooly <- c(as.vector(object$f$data[,2]),as.vector(object$g$data[,2]))
	datapooly <- rep(datapooly,c(object$f$counts,object$g$counts))
	
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


		if(is.null(h0sim)){	
			h0f <- object$f$globalH
			h0g <- object$g$globalH
			
		} else if(is.function(h0sim)){
			h0 <- h0sim(datapool[casetemp,],datapool[contemp,])

			if(!is.numeric(h0)) stop("if a function, 'h0sim' must return a numeric vector of length 1 or 2")
			
			if(length(h0)==1){
				h0f <- h0g <- h0
			} else if(length(h0)==2){
				h0f <- h0[1]
				h0g <- h0[2]
			} else {
				stop("if a function, 'h0sim' must return a numeric vector of length 1 or 2")
			}
						
		} else if(is.numeric(h0sim)){
			if(length(h0sim)!=2 && length(h0sim)!=1) stop("if numeric, 'h0sim' must be a vector of length 1 or 2")

			if(length(h0sim)==1){
				h0f <- h0g <- h0sim
			} else {
				h0f <- h0sim[1]
				h0g <- h0sim[2]
			}

		} else {
			stop("invalid 'h0sim' argument")
		}
		
		
		casep <- bivariate.density(data=datapool[casetemp,],pilotH=hpf,adaptive=F,res=res,comment=F,WIN=WIN)
		conp <- bivariate.density(data=datapool[contemp,],pilotH=hpg,adaptive=F,res=res,comment=F,WIN=WIN)
		
		casedens <- bivariate.density(data=datapool[casetemp,],pdef=casep,globalH=h0f,adaptive=T,res=res,WIN=WIN,gamma=exp((1/(n1+n2))*sum(c(log(1/sqrt(casep$zSpec)),log(1/sqrt(conp$zSpec))))),comment=F)
		condens <- bivariate.density(data=datapool[contemp,],pdef=conp,globalH=h0g,adaptive=T,res=res,WIN=WIN,gamma=exp((1/(n1+n2))*sum(c(log(1/sqrt(casep$zSpec)),log(1/sqrt(conp$zSpec))))),comment=F)
		
		risktemp <-risk(casedens,condens,log=object$log,plotit=F)$rsM
			
		mcmat_upper <- mcmat_upper+(risktemp>=rmat)
		if(comment) setTxtProgressBar(pb,i)
	}
	if(comment) close(pb)
	
	return(mcmat_upper/ITER)
}
