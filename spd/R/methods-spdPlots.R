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
#################################################################################
# S4 plot method
.plot.GPDTAILS<-function(x, which="ask",...)
{
if(missing(which)) which="ask"
	.interactivegpdtailsPlot(x, choices = c("Semi Parametric CDF", "Semi Parametric PDF",
					"Fit Assessment-Upper Tail", "Fit Assessment-Lower Tail"),
			plotFUN = paste(".plot.gpdtails.", 1:4, sep = ""), which = which, ...)
	# Return Value:
	invisible(x)
}

setMethod("plot", signature(x="GPDTAILS",y="missing"), .plot.GPDTAILS)

setMethod("plot", signature(x="GPDFIT", y="missing"), .plot.GPDFIT)
#---------------------------------------------------------------------
.interactivegpdtailsPlot = function(object, choices, plotFUN, which, ...)
{
# Some cecks:
    if (length(choices) != length(plotFUN))
    stop("Arguments choices and plotFUN must be of same length.")
    if (length(which) > length(choices))
    stop("Arguments which has incorrect length.")
    if (length(which) > length(plotFUN))
    stop("Arguments which has incorrect length.")
    # Plot:
    if (is.numeric(which)) {
        Which = rep(FALSE, times = length(choices))
        Which[which] = TRUE
    }

    if (which[1] == "all") {
        Which = rep(TRUE, times = length(choices))
    }

    if (which[1] == "ask") {
        .multgpdtailsPlot(object, choices, plotFUN, ...)
    } else {
        for ( i in 1:length(choices) ) {
            FUN = match.fun(plotFUN[i])
            if (Which[i]) FUN(object)
        }
    }

    # Return Value:
    invisible(object)
}

.multgpdtailsPlot = function(object, choices, ...)
{
    pick = 1
    while (pick > 0)
	{
        pick = menu(
        ### choices = paste("plot:", choices)
        choices = paste("plot:",choices,sep=""),
        title = "\nMake a plot selection (or 0 to exit):")
        switch (pick,
        .plot.gpdtails.1(object,...),  .plot.gpdtails.2(object,...),  .plot.gpdtails.3(object,...),
		.plot.gpdtails.4(object,...))
    	}
}

.plot.gpdtails.1 <- function(object, ...)
{
    lower = object@ptails$lower
    upper = object@ptails$upper
    origData = sort(object@data)

    minProbability = pspd(min(origData),object)
    maxProbability = pspd(max(origData),object)

    par.orig <- par(no.readonly = TRUE)
	N=length(origData)
    pLowerTail = seq(minProbability, lower, length.out=200)
    # sample lower tail
    pUpperTail = seq(upper, maxProbability, length.out=200)
    # sample upper tail
    pInterior  = seq(lower, upper, length.out=200)
    # sample interior
	nOrg=as.integer(seq(1,N,length.out=200))
    qLower = qspd(pLowerTail,object)
    qUpper = qspd(pUpperTail,object)
    qInterior = qspd(pInterior,object)
    xl=c(min(qLower),max(qUpper))
    yl=c(min(pLowerTail),max(pUpperTail))
	plot(sort(origData)[nOrg],ppoints(sort(origData)[nOrg]),type="p", col="steelblue",ylab="Probability",xlab="Returns",main="Semi-Parametric CDF")
    lines(qLower,pLowerTail,lwd=2,xlim=xl,ylim=yl, col="red")
    lines(qInterior,pInterior, lwd=3,col="darkgrey")
    lines(qUpper,pUpperTail, lwd=2,col="blue")
    leg.txt<-c("Pareto Lower Tail","Kernel Smoothed Interior","Pareto Upper Tail","Observed Data")
    legend(x="topleft", yjust=0,legend=leg.txt,lwd=3, lty=c(1,1,1,3), col=c('red','darkgrey','blue','steelblue'), cex=0.8)
	mtext(paste("spd  : GPD Tail Fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.7)
    par<-par.orig
    grid()
	box()
}


.plot.gpdtails.2 <- function(object, ...)
{
	lowerX = object@threshold$lower
	upperX = object@threshold$upper
	origData = sort(object@data)
	minX = min(origData)-0.02
	maxX = max(origData)+0.02
	
	par.orig <- par(no.readonly = TRUE)
	N=length(origData)
	dLowerTail = seq(minX, lowerX, length.out=200)
	# sample lower tail
	dUpperTail = seq(upperX, maxX, length.out=200)
	# sample upper tail
	dInteriorX  = seq(lowerX, upperX, length.out=200)
	# sample interior
	edf=density(origData,n=600)
	dLower = dspd(dLowerTail,object)
	dUpper = dspd(dUpperTail,object)
	dInterior = dspd(dInteriorX,object)
	xl=c(min(origData),max(origData))
	# max(dInterior) for certain contains max...usually, but we make sure
	yl=c(0,max(c(dLower,dInterior,dUpper)))
	plot(edf$x,edf$y,type="p", col="steelblue",ylab="Frequency",xlab="x",main="Semi-Parametric PDF")
	lines(dLowerTail,dLower,type="l",lwd=2,xlim=xl,ylim=yl, col="red")
	lines(dInteriorX,dInterior,lwd=3,col="darkgrey")
	lines(dUpperTail,dUpper, lwd=2,col="blue")
	leg.txt<-c("Pareto Lower Tail","Kernel Smoothed Interior","Pareto Upper Tail","Observed Data (Kernel Estimate)")
	legend(x="topleft", yjust=0,legend=leg.txt,lwd=3, lty=c(1,1,1,3), col=c('red','darkgrey','blue','steelblue'), cex=0.8)
	mtext(paste("spd  : GPD Tail Fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.7)
	par<-par.orig
	grid()
	box()
}

.plot.gpdtails.3 <- function(object, ...)
{
    par.orig <- par(no.readonly = TRUE)
	par(mfrow=c(2,2),oma=c(2,0,2,0))
    plot(object@fit$upperFit,which="all")
	title("GPD Fit Assessment-Upper Tail",outer=TRUE,cex=0.85)
    par<-par.orig
}

.plot.gpdtails.4<- function(object, ...)
{
    par.orig <- par(no.readonly = TRUE)
	par(mfrow=c(2,2),oma=c(2,0,2,0))
    plot(object@fit$lowerFit,which="all")
	title("GPD Fit Assessment-Lower Tail",outer=TRUE,cex=0.85)
    par<-par.orig
}