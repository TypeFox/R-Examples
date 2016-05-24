rsmc.sym <- function(object,symchoice,hpsim,h0sim,ITER,corrgrid,comment){
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
		
		if(symchoice=="f"){
			dat <- casetemp
		} else if(symchoice=="g"){
			dat <- contemp
		} else {
			dat <- range.n
		}
		
		if(is.null(hpsim)){
			hp <- object$f$pdef$pilotH
		} else if(is.function(hpsim)){
			hp <- hpsim(datapool[casetemp,],datapool[contemp,])
			
			if(!is.numeric(hp)) stop("if a function, 'hpsim' must return a numeric vector of length 1 or 2")
			
			if(length(hp)!=1&&length(hp)!=2){
				stop("if a function, 'h0sim' must return a numeric vector of length 1 or 2")
			} else if (length(hp)==2){
				hp <- hp[1]
				warning("'hpsim' length == 2, using first value only for symmetric pilot density bandwidth")
			}

		} else if(is.numeric(hpsim)){
			hp <- hpsim
			if(length(hp)>1){
				hp <- hp[1]
				warning("'hpsim' length > 1, using first value only for symmetric pilot density bandwidth")
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
		
		if(length(dat)==(n1+n2)){
			pil <- bivariate.density(data=datapool[dat,],pilotH=hp,adaptive=F,res=res,WIN=WIN,comment=F)
			gam <- exp(mean(log(1/sqrt(c(pil$zSpec)))))
		} else {
			pil <- bivariate.density(data=datapool[dat,],pilotH=hp,adaptive=F,res=res,WIN=WIN,comment=F,atExtraCoords=datapool[-dat,])
			gam <- exp(mean(log(1/sqrt(c(pil$zSpec,pil$zExtra)))))
		}
		
		casedens <- bivariate.density(data=datapool[casetemp,],pdef=pil,globalH=h0f,adaptive=T,res=res,WIN=WIN,gamma=gam,comment=F)
		condens <- bivariate.density(data=datapool[contemp,],pdef=pil,globalH=h0g,adaptive=T,res=res,WIN=WIN,gamma=gam,comment=F)
		
		risktemp <-risk(casedens,condens,log=object$log,plotit=F)$rsM
			
		mcmat_upper <- mcmat_upper+(risktemp>=rmat)
		if(comment) setTxtProgressBar(pb,i)
	}
	if(comment) close(pb)

	return(mcmat_upper/ITER)
}
