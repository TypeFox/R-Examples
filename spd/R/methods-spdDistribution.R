#################################################################################
##
##   R package spd by Alexios Ghalanos Copyright (C) 2008-2013
##   This file is part of the R package spd.
##
##   The R package spd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package spd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##	 Part of the psdp, qspd functions inspired by Rene Carmona's book
##   "STATISTICAL ANALYSIS of FINANCIAL DATA in S-PLUS"
#################################################################################

# Developer Notes: p,d,q methods switch to sub methods given tail distribution (e.g.gpd)
# via the fit class (first 3 letter) thereby enabling easy extension of methods to other
# tail distributions.
# ------------------------------------------------------------------------------
pspd = function(q, fit, linear=TRUE)
{
	UseMethod("pspd")
}

.pspd = function(q, fit, linear = TRUE)
{
	tail=substr(is(fit)[1],1,3)
	ans<-switch(tail,
			GPD=.pspdgpd(q, fit, linear))
	return(ans)
}
.pspdgpd = function(q, fit, linear = TRUE)
{
    x <- sort(q)
    N <- length(x)
    val <- vector(length = N,mode = "numeric")
    origData <-sort(as.matrix(fit@data))
    n <- length(origData)

    # define regions of parametric and kernel interaction
    lowerThresh = fit@threshold$lower
    upperThresh = fit@threshold$upper
    upperPar = fit@fit$upperFit@fit$par.est
    lowerPar = fit@fit$lowerFit@fit$par.est
    # interior coordinates

    # this is the estimate of CDF in the interval
    pp = fit@fit$kernelFit
    ppx = pp$x
    ppy = pp$y
    ppz = cumsum(pp$y/sum(pp$y))
	
	if(any(x>max(ppx)) && any(x<min(ppx))) {
		kf = bkde(origData, fit@kernel, gridsize=as.integer(length(origData)),range.x=c(min(x),max(x)))
		ppx = kf$x
		ppy = kf$y
		ppz = cumsum(kf$y/sum(kf$y))
	}
	if(any(x>max(ppx)) && !any(x<min(ppx))) {
		kf = bkde(origData, fit@kernel, gridsize=as.integer(length(origData)),range.x=c(min(ppx),max(x)))
		ppx = kf$x
		ppy = kf$y
		ppz = cumsum(kf$y/sum(kf$y))
	}
	if(!any(x>max(ppx)) && any(x<min(ppx))) {
		kf = bkde(origData, fit@kernel, gridsize=as.integer(length(origData)),range.x=c(min(x),max(ppx)))
		ppx = kf$x
		ppy = kf$y
		ppz = cumsum(kf$y/sum(kf$y))
	}
	
	# interior values
	midx <-  lowerThresh <= x & x <= upperThresh
    xm <- as.numeric(x[midx])
    nm <- length(xm)

	# The kernel estimation is done on all the data points so we must truncate...
    # This is the count of the number of kernel fitted occurences < Lower Threshold
    midstart <- sum(ppx < lowerThresh)

    # This is the count of the number of kernel fitted occurences <= Upper Threshold
    midend <- sum(ppx <= upperThresh) + 1

    # This is the kenrel fitted data in interior
    midpts <- as.numeric(ppx[midstart:midend])

    # This is the number of kernel fitted points in between
	nmidpts <- length(midpts)
    oldind <- vector(length = nm,mode = "numeric")
    oldind[1:nm] <- -1

    # Now we find the location of xm (midpoints of p input data) in the midpoints region of the
    # kernel fitted data
    options(show.error.messages = FALSE)
    midindex <- approx(x=sort(midpts),y=1:nmidpts, xout=xm,method="linear",ties="ordered")$y+ midstart
    options(show.error.messages = TRUE)

    # ppx data corresponding to the p data
    midval <- ppx[midindex]
	midval[is.na(midval)] <- 0
	midindexL <- midindex -1
	midindexL[midindexL == 0] <- NA
    # shifted [backwards] data corresponding to the p data for linear interpolation
    midvalS <- ppx[midindexL]
	midvalS[is.na(midvalS)] <- (x[midx])[is.na(midvalS)]
    # corresponding cumulative probability ppz
    pmidval <- ppz[midindex]
	pmidval[is.na(pmidval)] <- 0

	pmidvalS <- ppz[midindexL]
	pmidvalS[is.na(pmidvalS)] <- 0

    if (linear)
    {
    	val[midx]<- pmidvalS + (pmidval - pmidvalS)*(xm-midvalS)/(midval - midvalS)
    }
    else
    {
    	val[midx] <- pmidvalS
    }
    # this is the estimate of CDF at upper threshold:
    p.upper <- ppz[midend - 1] + (ppz[midend]-ppz[midend - 1] )*(upperThresh - ppx[midend-1])/(ppx[midend]- ppx[midend-1])
    p.lower <- ppz[midstart] + (ppz[midstart + 1]-ppz[midstart] )*(lowerThresh - ppx[midstart])/(ppx[midstart + 1]- ppx[midstart])
    upper.tail.x <- x > upperThresh
    lower.tail.x <- x < lowerThresh
    k <- upperPar[1]
    a <- upperPar[2]
    b <- upperThresh
    if(length(x[x>upperThresh])!=0)
    {
    	val[upper.tail.x] <- 1-(1-p.upper)*pgpd(x[x>upperThresh]-upperThresh,xi=upperPar[1],mu=0,beta=upperPar[2],lower.tail=FALSE)
	    if (k < 0 & sum(x > b - a/k) > 0)
	    {
	    	val[x > b - a/k]  <- 1.0
	    }
    }
    if(length(x[x<lowerThresh])!=0)
    {
    # this is the estimate of CDF at lower threshold:
	    k <- lowerPar[1]
	    a <- lowerPar[2]
	    b <- lowerThresh
		if(length(x[x<lowerThresh])==0) xlin=0 else xlin=x[x<lowerThresh]
	    val[lower.tail.x] <- p.lower*pgpd(lowerThresh-xlin,xi=lowerPar[1],mu=0,beta=lowerPar[2],lower.tail=FALSE)
	    if (k < 0 & sum(x < b + a/k) > 0)
	    {
	    	val[x < b + a/k]<-0
	    }
    }
    #val <- approx(x=ppx,y=val, xout=sort(q),method="linear",ties="ordered")$y
	retval<-val
    retval[sort.list(q)] <- val
    return(retval)
}


setMethod("pspd", signature(q="vector", fit="SPD", linear="logical"), .pspd)
setMethod("pspd", signature(q="vector", fit="SPD", linear="missing"), .pspd)
# ------------------------------------------------------------------------------
dspd <- function(x, fit, linear=TRUE)
{
	UseMethod("dspd")
}

.dspd = function(x, fit, linear = TRUE)
{
	tail=substr(is(fit)[1],1,3)
	ans<-switch(tail, GPD=.dspdgpd(x, fit, linear))
	return(ans)
}

.dspdgpd <- function(x, fit, linear = TRUE)
{
    X <- sort(x)
    N <- length(x)
	tpm=vector(mode="numeric")
    val <- vector(length = N,mode = "numeric")
    origData <-sort(as.matrix(fit@data))
    n <- length(origData)
    lowerThresh = fit@threshold$lower
    upperThresh = fit@threshold$upper
    upperPar = fit@fit$upperFit@fit$par.est
    lowerPar = fit@fit$lowerFit@fit$par.est

    # this is the estimate of CDF in the interval
    pp = fit@fit$kernelFit
    ppx = pp$x
    ppy = pp$y
    ppz = cumsum(pp$y/sum(pp$y))
	if(any(x>max(ppx)) && any(x<min(ppx))) {
		kf = bkde(origData, fit@kernel, gridsize=as.integer(length(origData)),range.x=c(min(x),max(x)))
		ppx = kf$x
		ppy = kf$y
		ppz = cumsum(kf$y/sum(kf$y))
	}
	if(any(x>max(ppx)) && !any(x<min(ppx))) {
		kf = bkde(origData, fit@kernel, gridsize=as.integer(length(origData)),range.x=c(min(ppx),max(x)))
		ppx = kf$x
		ppy = kf$y
		ppz = cumsum(kf$y/sum(kf$y))
	}
	if(!any(x>max(ppx)) && any(x<min(ppx))) {
		kf = bkde(origData, fit@kernel, gridsize=as.integer(length(origData)),range.x=c(min(x),max(ppx)))
		ppx = kf$x
		ppy = kf$y
		ppz = cumsum(kf$y/sum(kf$y))
	}

    midx <-  lowerThresh <= ppx & ppx <= upperThresh
    xm <- as.numeric(ppx[midx])
    nm <- length(xm)
    # This is the count of the number of kernel fitted occurences < Lower Threshold
	midstart <- sum(ppx < lowerThresh)

    # This is the count of the number of kernel fitted occurences <= Upper Threshold
	midend <- sum(ppx <= upperThresh) + 1
    # this is the estimate of PDF in the interval:
    valtry<-try(xg<-.interpPDF(xm, ppx, ppz), silent=TRUE)
	if (!inherits(valtry,"try-error"))
	{
		val[midx]<-.interpPDF(xm, ppx, ppz)
	}
    # this is the estimate of F at upper threshold:
    upper.tail.x <- ppx > upperThresh
    lower.tail.x <- ppx < lowerThresh
    if(length(ppx[ppx>upperThresh])!=0)
    {
	    k <- upperPar[1]
	    a <- upperPar[2]
	    b <- upperThresh
		start=val[midx][length(val[midx])]
		#func = function(x) (1-p.upper)*pgpd(ppx,xi=k,mu=0,beta=a)
		zz=dgpd(ppx[ppx>upperThresh]-upperThresh,xi=k,mu=0,beta=a)
		xd=function(x) (zz/sum(zz)*((1-x)*sum(ppy)))[1]-start
		ff=uniroot(f=xd,interval=c(0.5,1))
	    val[upper.tail.x] <- zz/sum(zz)*((1-ff$root)*sum(ppy))
    }
    if(length(ppx[ppx<lowerThresh])!=0)
    {
	    # this is the estimate of F at lower threshold:
	    k <- lowerPar[1]
	    a <- lowerPar[2]
	    b <- lowerThresh
		zz=dgpd((lowerThresh-ppx[ppx<lowerThresh]),xi=k,mu=0,beta=a)
		start=val[midx][1]
		xd=function(x) (zz/sum(zz)*(x*sum(ppy)))[length(zz)]-start
		ff=uniroot(f=xd,interval=c(0,0.5))
	    val[lower.tail.x] <- zz/sum(zz)*(ff$root*sum(ppy))
    }
	val <- approx(x=ppx,y=val, xout=sort(x),method="linear",ties="ordered")$y
	retval<-val
    retval[sort.list(x)] <- val
    return(retval)
}


setMethod("dspd", signature(x="vector", fit="SPD", linear="logical"), .dspd)
setMethod("dspd", signature(x="vector", fit="SPD", linear="missing"), .dspd)

# ------------------------------------------------------------------------------
qspd = function(p, fit, linear=TRUE)
{
	UseMethod("qspd")
}

.qspd = function(p, fit, linear = TRUE)
{
	tail=substr(is(fit)[1],1,3)
	ans<-switch(tail,
			GPD=.qspdgpd(p, fit, linear))
	return(ans)
}


.qspdgpd<-function(p,fit,linear = TRUE)
{
    x = sort(p)
    N <- length(x)
    origData <-sort(as.matrix(fit@data))
    n <- length(origData)
    # here we need to check: IF GDP ok, else need to calculate upper/lower Thresh
    lowerThresh = fit@threshold$lower
    upperThresh = fit@threshold$upper
    upperPar = fit@fit$upperFit@fit$par.est
    lowerPar = fit@fit$lowerFit@fit$par.est
    # this is the estimate of F in the interval:
    val <- vector(length = N,mode = "numeric")
    validp <- x >=0 & x <= 1
    val[!validp] <- NA
	midstart <- sum(origData < lowerThresh)
	midend <- sum(origData <= upperThresh) + 1
    pp <- ppoints(origData)
    p.upper <- pp[midend - 1] + (pp[midend]-pp[midend - 1] )*(upperThresh - origData[midend-1])/(origData[midend]- origData[midend-1])
    p.lower <- pp[midstart] + (pp[midstart + 1]-pp[midstart] )*(lowerThresh - origData[midstart])/(origData[midstart + 1]- origData[midstart])
	midx <-  p.lower <= x & x <= p.upper
    xm <- as.double(x[midx])
    nm <- as.integer(length(xm))
    midpts <- as.numeric(pp[midstart:midend])
	nmidpts <- length(midpts)
    oldind <- vector(length = nm,mode = "numeric")
    oldind[1:nm] <- -1
	midindex <- approx(y=1:nmidpts, x=midpts, xout=xm, method="constant", ties="ordered")$y + midstart
	midval <- pp[midindex]

	midindexL <- midindex -1
	midindex[midindex == 0] <- NA
	midindex[is.na(midindex)] <- 0

	midvalS <- pp[midindexL]
	midvalS[is.na(midvalS)] <- 0

	pmidval <- origData[midindex]
	pmidval[is.na(pmidval)] <- 0

	pmidvalS <- origData[midindexL]
	pmidvalS[is.na(pmidvalS)] <- 0

    if (linear)
    {
    	val[midx] <- pmidvalS + (pmidval - pmidvalS)*(xm-midvalS)/(midval - midvalS)
    }
    else
    {
    	val[midx] <- pmidvalS
    }
    # this is the estimate of F at upper threshold:
    upper.tail.x <- x > p.upper & validp
    lower.tail.x <- x < p.lower & validp
    if(length(val[upper.tail.x])!=0)
    {
	    xu=x[x > p.upper & validp]
	    k <- upperPar[1]
	    a <- upperPar[2]
	    b <- upperThresh
	    val[upper.tail.x] <- b + qgpd((xu-p.upper)/(1-p.upper), xi=k, mu=0, beta=a)
    }
    # this is the estimate of F at lower threshold:
    if(length(val[lower.tail.x])!=0)
    {
	    xl = x[x < p.lower & validp]
	    k <- lowerPar[1]
	    a <- lowerPar[2]
	    b <- lowerThresh
	    val[lower.tail.x] <- b - qgpd(1-xl/p.lower, xi=k, mu=0, beta=a)
    }
    retval <- val
    retval[sort.list(p)] <- val
    return(retval)
}

setMethod("qspd", signature(p="vector", fit="SPD",linear="logical"), .qspd)
setMethod("qspd", signature(p="vector", fit="SPD",linear="missing"), .qspd)
# ------------------------------------------------------------------------------
rspd = function(n, fit, linear=TRUE)
{
	UseMethod("rspd")
}

.rspd<-function(n, fit, linear=TRUE)
{
	# if(rseed!=0) set.seed(rseed)
    # use the inverse CDF method
    z<-runif(n)
    retval<-.qspd(z, fit, linear)
    return(retval)
}
setMethod("rspd", signature(n="vector", fit="SPD",linear="logical"), .rspd)
setMethod("rspd", signature(n="vector", fit="SPD",linear="missing"), .rspd)

# ------------------------------------------------------------------------------



