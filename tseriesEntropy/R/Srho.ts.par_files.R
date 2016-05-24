
  ## tseriesEntropy
  ## Entropy based tests of serial dependence and nonlinearity
  #
  #  The authors of this software is
  #
  #  Simone Giannerini, Copyright (c) 2009.
  #
  #  Permission to use, copy, modify, and distribute this software for any
  #  purpose without fee is hereby granted, provided that this entire notice
  #  is included in all copies of any software which is or includes a copy
  #  or modification of this software and in all copies of the supporting
  #  documentation for such software.
  #
  #  This program is free software; you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation; either version 2 of the License, or
  #  (at your option) any later version.
  #
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.
  #
  #  A copy of the GNU General Public License is available at
  #  http://www.r-project.org/Licenses/

## **************************************************************************************************
## PARALLEL VERSION OF THE TESTS
## **************************************************************************************************

Srho.test.AR.p <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference","mlcv", "lscv"), method =c("integral","summation"), maxpts=0, tol=1e-03,
order.max=10, fit.method=c("yule-walker", "burg", "ols", "mle", "yw"),smoothed=TRUE,
nslaves=detectCores()){
	if(nslaves < 2) stop('Number of slaves must be at least 2.')
	if(nslaves > detectCores()) warning('Number of slaves greater than number of cores.')
	
	if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
	if(length(quant)==1){
		if(quant==0.99){
		quant <- c(0.95,quant)
		}else{
		quant <- c(quant,0.99)
		}
	}
    Bnew       <- B + 2*nslaves # adds 2 replications per slave
    bw         <- match.arg(bw)
	method     <- match.arg(method)
	fit.method <- match.arg(fit.method)
    arg.s      <- list(x=x,order.max=order.max,fit.method=fit.method,nsurr = Bnew)
    fun        <- switch(smoothed+1,"surrogate.AR","surrogate.ARs")
	if (missing(y)){
        x.surr <- do.call(eval(fun),args=arg.s)
        cl     <- makeCluster(nslaves)
        junk1  <- clusterEvalQ(cl, library(tseriesEntropy))
        varlist <- c("x.surr","B","lag.max","bw","method","fit.method","maxpts","tol","arg.s")
        for(v in varlist) {
            environment(v) <- .GlobalEnv
        }
        clusterExport(cl,varlist,envir = environment())
        dum <- clusterSplit(cl,1:Bnew)
        result <- clusterApply(cl,x=dum,fun=function(x){
            safe.Srho(x.surr=x.surr$surr[,x],B.good=(length(x)-2),lag.max=lag.max,bw=bw,method=method, maxpts=maxpts, tol=tol)
            })
#		result   <-  clusterEvalQ(cl,tseriesEntropy:::safe.Srho(x.surr=x.surr$surr,B.good=B,lag.max=lag.max,bw=bw,method=method, maxpts=maxpts, tol=tol));
		stopCluster(cl);
		S.x    <- Srho.ts(x,lag.max=lag.max,bw=bw,method=method, plot=FALSE,
					maxpts=maxpts, tol=tol)
		M <- NULL
		for(i in 1:nslaves){
			M <- cbind(M,result[[i]])
		}
        rownames(M) <- 1:lag.max
        M.95 <- apply(M,1,quantile,probs=c(quant[1]));
		M.99 <- apply(M,1,quantile,probs=c(quant[2]));
#        names(S.x) <- 1:lag.max
        ind95  <- which(S.x>=M.95);
        ind99  <- which(S.x>=M.99);
        out            <- new("Srho.test")
        out@.Data      <- S.x@.Data
        out@lags       <- 1:lag.max
        out@stationary <- FALSE
        out@data.type  <- "continuous"
        out@test.type  <- "nonlinearity (AR)"
        out@call       <- match.call();
        out@call.h     <- x.surr$call
        out@quantiles  <- cbind(M.95,M.99)
        q.names        <- paste("Q",as.character((quant*100)),"%",sep='')
        colnames(out@quantiles)     <- q.names
        rownames(out@quantiles)     <- 1:lag.max
        out@significant.lags        <- list(as.integer(names(ind95)),as.integer(names(ind99)))
        names(out@significant.lags) <- q.names
        out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
        names(out@p.value)          <- 1:lag.max
        if (plot) {
            plot(out)
            return(invisible(out))
        }
        else return(out)

    } else {
        return(cat("Cross-Entropy testing for non linearity not yet implemented"))
    }
}
## **************************************************************************************************

Trho.test.SA.p <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference","mlcv", "lscv"), method =c("integral","summation"), maxpts=0, tol=1e-03,
nlag=trunc(length(x)/4),Te=0.0015,RT=0.9,eps.SA=0.01,nsuccmax=30,nmax=300,che=100000,
nslaves=detectCores()){
	if(nslaves < 2) stop('Number of slaves must be at least 2.')
	if(nslaves > detectCores()) warning('Number of slaves greater than number of cores.')
	
	if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
	if(length(quant)==1){
		if(quant==0.99){
		quant <- c(0.95,quant)
		}else{
		quant <- c(quant,0.99)
		}
	}
    Bnew       <- B + 2*nslaves # adds 2 replications per slave
    bw         <- match.arg(bw)
	method     <- match.arg(method)
    arg.s      <- list(x=x,nlag=nlag,nsurr=Bnew,Te=Te,RT=RT,eps.SA=eps.SA,nsuccmax=nsuccmax,nmax=nmax,che=che)
    fun        <- "surrogate.SA"
	if (missing(y)){
        x.surr <- do.call(eval(fun),args=arg.s)
        cl     <- makeCluster(nslaves)
        junk1  <- clusterEvalQ(cl, library(tseriesEntropy))
        varlist <- c("x.surr","B","lag.max","bw","che","maxpts","tol")
        for(v in varlist) {
            environment(v) <- .GlobalEnv
        }
        clusterExport(cl,varlist,envir = environment())
        dum <- clusterSplit(cl,1:Bnew)
        result <- clusterApply(cl,x=dum,fun=function(x){
            safe.Trho(x.surr=x.surr$surr[,x],B.good=(length(x)-2),lag.max=lag.max,bw=bw,method=method,maxpts=maxpts, tol=tol)
            })
#		result   <-  clusterEvalQ(cl,tseriesEntropy:::safe.Srho(x.surr=x.surr$surr,B.good=B,lag.max=lag.max,bw=bw,method=method, maxpts=maxpts, tol=tol));
		stopCluster(cl);
        S.x    <- (Srho.ts(x,lag.max=lag.max,bw=bw,method=method, plot=FALSE,maxpts=maxpts, tol=tol)@.Data
                         - Srho.cor(x,lag.max=lag.max,plot=FALSE)@.Data)^2
		M <- NULL
		for(i in 1:nslaves){
			M <- cbind(M,result[[i]])
		}
        rownames(M) <- 1:lag.max
		M.95 <- apply(M,1,quantile,probs=c(quant[1]));
		M.99 <- apply(M,1,quantile,probs=c(quant[2]));
#		names(S.x) <- 1:lag.max
		ind95 <- which(S.x>=M.95);
		ind99 <- which(S.x>=M.99);
		out            <- new("Srho.test")
		out@.Data      <- S.x@.Data
		out@lags       <- 1:lag.max
		out@stationary <- FALSE
		out@data.type  <- "continuous"
		out@test.type  <- "nonlinearity (SA)"
		out@call       <- match.call();
#		out@call.h     <- x.surr$call;
		out@quantiles  <- cbind(M.95,M.99)
		q.names <- paste("Q",as.character((quant*100)),"%",sep='')
		colnames(out@quantiles)     <- q.names
		rownames(out@quantiles)     <- 1:lag.max
		out@significant.lags        <- list(as.integer(names(ind95)),as.integer(names(ind99)))
		names(out@significant.lags) <- q.names
        out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
        names(out@p.value)          <- 1:lag.max
		if (plot) {
			plot(out,ylab = "T")
			return(invisible(out))
		}
		else return(out)

    } else {
        return(cat("Cross-Entropy testing for non linearity not yet implemented"))
    }
}

## ****************************************************************************************************

Trho.test.AR.p <- function(x, y, lag.max = 10,  B = 100, plot = TRUE, quant = c(0.95, 0.99),
bw = c("reference","mlcv", "lscv"), method =c("integral","summation"), maxpts=0, tol=1e-03,
order.max=10, fit.method=c("yule-walker", "burg", "ols", "mle", "yw"), smoothed=TRUE,
nslaves=detectCores()){
    if(nslaves < 2) stop('Number of slaves must be at least 2.')
	if(nslaves > detectCores()) warning('Number of slaves greater than number of cores.')

	if(any(quant<=0|quant>=1)) stop("elements of quant must lie in ]0,1[");
	if(length(quant)==1){
		if(quant==0.99){
		quant <- c(0.95,quant)
		}else{
		quant <- c(quant,0.99)
		}
	}
    Bnew       <- B + 2*nslaves # adds 2 replications per slave
    bw         <- match.arg(bw)
	method     <- match.arg(method)
	fit.method <- match.arg(fit.method)
    arg.s      <- list(x=x,order.max=order.max,fit.method=fit.method,nsurr = Bnew)
    fun        <- switch(smoothed+1,"surrogate.AR","surrogate.ARs")
	if (missing(y)){
        x.surr <- do.call(eval(fun),args=arg.s)
        cl     <- makeCluster(nslaves)
        junk1  <- clusterEvalQ(cl, library(tseriesEntropy))
        varlist <- c("x.surr","B","lag.max","bw","method","fit.method","maxpts","tol","arg.s")
        for(v in varlist) {
            environment(v) <- .GlobalEnv
        }
        clusterExport(cl,varlist,envir = environment())
        dum <- clusterSplit(cl,1:Bnew)
        result <- clusterApply(cl,x=dum,fun=function(x){
            safe.Trho(x.surr=x.surr$surr[,x],B.good=(length(x)-2),lag.max=lag.max,bw=bw,method=method,maxpts=maxpts, tol=tol)
            })
#		result   <-  clusterEvalQ(cl,tseriesEntropy:::safe.Srho(x.surr=x.surr$surr,B.good=B,lag.max=lag.max,bw=bw,method=method, maxpts=maxpts, tol=tol));
		stopCluster(cl);

        S.x    <- (Srho.ts(x,lag.max=lag.max,bw=bw,method=method, plot=FALSE,
                                maxpts=maxpts, tol=tol)@.Data - Srho.cor(x,lag.max=lag.max,plot=FALSE)@.Data)^2
		M <- NULL
		for(i in 1:nslaves){
			M <- cbind(M,result[[i]])
		}
        rownames(M) <- 1:lag.max
		M.95 <- apply(M,1,quantile,probs=c(quant[1]));
		M.99 <- apply(M,1,quantile,probs=c(quant[2]));
#		names(S.x) <- 1:lag.max
		ind95 <- which(S.x>=M.95);
		ind99 <- which(S.x>=M.99);
		out            <- new("Srho.test")
		out@.Data      <- S.x
		out@lags       <- 1:lag.max
		out@stationary <- FALSE
		out@data.type  <- "continuous"
		out@test.type  <- "nonlinearity (AR)"
		out@call       <- match.call();
		#out@call.h     <- NULL
		out@quantiles  <- cbind(M.95,M.99)
		q.names <- paste("Q",as.character((quant*100)),"%",sep='')
		colnames(out@quantiles) <- q.names
		rownames(out@quantiles) <- 1:lag.max
		out@significant.lags    <- list(as.integer(names(ind95)),as.integer(names(ind99)))
		names(out@significant.lags) <- q.names
        out@p.value                 <- rowMeans(M >= S.x) # bootstrap p-value
        names(out@p.value)          <- 1:lag.max
		if (plot) {
		plot(out,ylab = "T")
		return(invisible(out))
		}
		else return(out)
	
	} else {
		return(cat("Cross-Entropy testing for non linearity not yet implemented"))
	}
}
## **************************************************************************************************
