##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2013 Karim Rahim 
##
##     Written by Karim Rahim and Wesley Burr.
##
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com


##################################################################
##
##  plot.mtm
##
##  Takes a mtm object, and plots either the associated spectrum
##  (obj$spec) or the harmonic F-test statistic (obj$Ftest).
##
##################################################################
plot.mtm <- function(x, 
                     jackknife=FALSE, 
                     Ftest=FALSE, 
                     ftbase=1.01,
                     siglines=NULL, 
                     ...) {
    
    ## Set frequency axis and label
    dtUnits <- x$mtm$dtUnits
    deltaT <- x$mtm$deltaT
    
    ## if the user has not set 'xlab' ... set it for them:
    if(!hasArg("xlab")) {
        if(!(x$mtm$dtUnits == "default")) {
            xlab <- paste("Frequency in cycles/",dtUnits,sep="") }
        else {
            xlab <- paste("Frequency")
        }
    } 
    
    if(Ftest) {
        if(!hasArg("xlab")) {
        .plotFtest(x,xlab=xlab,siglines=siglines,ftbase=ftbase, ...)
    } else {
        .plotFtest(x, siglines=siglines, ftbase=ftbase, ...)
    }  
    } 
    else 
    { ## plot spectrum only
        ## modified to remove calls to plot.spec
        ## for R version 3.1.0
        ##
        class(x) <- "spec"
        if(x$mtm$taper=="sine") {
            if(!hasArg("xlab")) {
                plot( x, xlab=xlab, sub=" ", ...)
            } else {
                plot( x, sub=" ", ...) 
            }  
        }
        else { ## case of taper=="dpss"
            nw <- x$mtm$nw
            k <- x$mtm$k
            sub <- paste("(NW = ", nw, " K = ", k,")", sep="")
            log <- match.call(expand.dots = )$log
            if(jackknife) {
                dBPlot <- FALSE
                if(!is.null(log) && log== "dB" ) {
                    dBPlot <- TRUE }
                if(jackknife && !is.null(x$mtm$jk)) {
                    if(dBPlot) {
                        upperCI <- 10*log10(x$mtm$jk$upperCI)
                        lowerCI <- 10*log10(x$mtm$jk$lowerCI)
                        minVal <- 10*log10(x$mtm$jk$minVal)
                        maxVal <- 10*log10(x$mtm$jk$maxVal) 
                    } 
                    else {
                        upperCI <- x$mtm$jk$upperCI
                        lowerCI <- x$mtm$jk$lowerCI
                        minVal <- x$mtm$jk$minVal
                        maxVal <- x$mtm$jk$maxVal
                    }
                    yRange <- c(minVal, maxVal)
                    if(!hasArg("xlab")) {
                        .lplotSpec( x, xlab=xlab, sub=sub, ylim=yRange, ...)
                    } else {
                        .lplotSpec( x, sub=sub, ylim=yRange, ...)
                    }  
                    lines(x$freq, upperCI, lty=2, col=2)
                    lines(x$freq, lowerCI, lty=2, col=3)
                }
            }
            else {
                if(!hasArg("xlab")) {
                    .lplotSpec( x, xlab=xlab, sub=sub, ...) 
                } else {
                    .lplotSpec( x, sub=sub, ...)
                }
            } 
        } ## end of dpss case
    } ## spectrum plot end
} ## end of function

##################################################################
##
##  plot.mtm.coh
##
##  Takes a mtm.coh object, and plots the Magnitude-Squared 
##  Coherence, with multiple y-axes.
##
##################################################################
plot.mtm.coh <- function(x, 
                         percentGreater=NULL,
                         nehlim=10, 
                         nehc=4,
                         cdfQuantilesTicks=NULL,
                         drawPercentLines=TRUE,
                         percentG=c(.1,.2,.5,.8,.9), 
                         ...) {

    if(  is.null(x$NTmsc) || is.null(x$NTvar)  || is.null(x$msc)
       || is.null(x$freq) || is.null(x$nfreqs) || is.null(x$k)) {
        stop("Requires mtm.coh object. Run mtm.coh on two mtm objects with returnInternals=TRUE.")
    }
    
    TRmsc <- x$NTmsc
    NTvar <- x$NTvar
    freqs <- x$freq
    nfreqs <- x$nfreqs
    k <- x$k
    
    ##nehlim and nehc are for smoothing 
    ## currently we plot the smoothed transformed coherence
    ## and lower CI after smoothing the variance
    plotTRmsc <- .lftr3p(TRmsc, NTvar, nfreqs,
                       nehlim,nehc, "even", "ext")
    trnrm_ <- .trnrm(k)
    par(oma=c(2,4,0,2))
    plot.new()
    ## note the ... was mainly implemented for xaxs="i"
    ## Undefined behaviour with other options 
    plot.window(range(freqs), range(plotTRmsc[,2]), ...)
    xy <- xy.coords(freqs,plotTRmsc[,2])
    ## plot smoothed msc
    plot.xy(xy, type="l", lwd=1, ...)
    ## plot one sd dev lower jackknife variance
    lines(freqs, plotTRmsc[,1], lty=3, lwd=1)
    box()
    axis(1)
    
    ## allow for user-settable xlabel, or unit display
    if(!hasArg("xlab")) {
        if(!(x$mtm$dtUnits == "default")) {
            xlab <- paste("Frequency in cycles/",x$mtm$dtUnits,sep="") }
        else {
            xlab <- paste("Frequency")
        }
    }
    mtext(xlab, side=1, line=3)
    
    ## basic left axis
    axis(2)
    mtext("Arctanh Transform of MSC",
          side=2, line=2, cex=par()$cex)
    
    ##  outer MSC axis on the left
    msc <- .FtoMSC(plotTRmsc[,2], trnrm_)
    mscTicks <- pretty(msc)
    
    
    ## transform ticks for at
    ##C2toF is coherence to inverse transform
    TRmscTicks <- .C2toF(mscTicks, trnrm_)
    axis(2, at=TRmscTicks, labels=mscTicks, outer=TRUE)
    mtext("Magnitude Squared Coherence", side=2, line=6)
    
    ##mscToCDF values may have issues for highly coherent values
    ## values over .9 will cause issues
    if(is.null(cdfQuantilesTicks)) {
        cdfQuantiles <- .mscToCDFquantiles(msc, k)
        cdfQuantilesTicks <- pretty(cdfQuantiles)
    }

    ## put right axis
    Qlvl <- .cdfToMSqCoh(cdfQuantilesTicks, k)
    TRQlvl <- .C2toF(Qlvl, trnrm_)
    
    cumulativeDistVals <- .C2toF(msc, trnrm_)
    axis(4, at=TRQlvl, labels=cdfQuantilesTicks)
    
    mtext("CDF for Independent Data",
          side=4, line=2) 

    if(drawPercentLines == TRUE) {
        percentGprob <- percentG
        percentG <- .C2toF(.cdfToMSqCoh(percentG, k),  trnrm_)
        lenPercentG <- length(percentG)
        for(i in 1:lenPercentG) {
            lines(freqs, array(percentG[i], nfreqs), lty=2)
        }
    }
    
    if(!is.null(percentGreater)) {
        mtext(paste("CDF for C=   10.0% 20.0% 50.0% 80.0% 90.0%"),
              side=1, line=4, adj=-1, cex=.8)
        mtext(paste("% of data > Q     ",
                    100*round( percentGreater[1], digits=3),
                    "% ",
                    100*round( percentGreater[2], digits=3),
                    "% ",
                    100*round( percentGreater[3], digits=3),
                    "% ",
                    100*round( percentGreater[4], digits=3),
                    "% ",
                    100*round( percentGreater[5], digits=3),
                    "%", sep=""),
              side=1, line=5, adj=-1, cex=0.8)
    }
    return(list(sigProb = percentGprob, sigNT = percentG))
}
