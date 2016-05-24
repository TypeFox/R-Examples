##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2013 Karim Rahim 
##
##     Written by Karim Rahim.
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
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com


## siglines degrees of freedom correction Oct 4, 2012 karim

##################################################################
##
##  .plotFtest
##
##  Takes a mtm object, and plots the harmonic F-test statistic 
##  (obj$Ftest).
##
##################################################################
.plotFtest <- function(x, 
                       ftbase=1.01, 
                       siglines=NULL, 
                       xlab="Frequency", 
                       ...) {

    if(is.null(x$mtm$Ftest) || !("Ftest" %in% class(x))) {
      stop(paste("Ftest not computed for given mtm object!"))
    }

    ## correct notes and warnings with devel R version to upload to cran.
    ##arglist <- list(...)
    ##log <- arglist$log ##match.call(expand.dots = )$log
    ##ylab <- arglist$ylab ##match.call(expand.dots = )$ylab
    log <- match.call(expand.dots = )$log
    ylab <- match.call(expand.dots = )$ylab
    
    if(is.null(ylab)) ylab <- "Harmonic F-test Statistic"
    
    ylog = "n"
    if(is.null(log) || log == "yes") {
        ylog = "y"
    }

    ftestVals = x$mtm$Ftest
    ftestVals[ftestVals < ftbase] <- ftbase
    ftmax <- max(ftestVals)

    .lplotDefault(x$freq, ftestVals, log=ylog, ylab=ylab, xlab=xlab, 
                 ylim=c(ftbase,ftmax), type="l", ...)
    
    ## add siglines if defined
    if(!is.null(siglines)) {
        for(j in 1:length(siglines)) {
            if(is.numeric(siglines[j]) && 0.80 <= siglines[j] && 1.000000 >= siglines[j]) {
                ## degree of freedom correction P&W page 499 changed to 2, 2*k-2 date Sept 30 2012
                sig0 <- qf(siglines[j],2,2*x$mtm$k-2) 
                abline(h=sig0, col="red", lty=2, lwd=1) 
                mtext(paste(siglines[j]*100,"%",sep=""), side=4, line=0, at=ceiling(sig0), col="red")
            }
        } ## end for
    } ## end logical
}

## this is a hack method to strip depreciated hidden parameters.
## This suggestion is from Gavin Simpson's blog and he traces it
## to a suggestion from Brian Ripley
## see: http://ucfagls.wordpress.com/2011/07/23/
## local hidden plotting routine to strip depreciated parameter
## dT from arguments lost
.lplotSpec <- function(x, ..., dT) {
    ## should call plot.spec prior to 3.1 dev.
    plot(x, ...)
}

## this is currently used in F--test plots.
## modified Feb 2015 to compile for cran
.lplotDefault <- function(x, y, ..., dT) {
    ## should call plot default 
    plot(x, y, ...)
}

## utilities functions added for plot.mtm.coh
## fortran versions exist and could be cleaned up and implemented...
## nhi is nfreqs
## c
## c x       Input Data
## c xv      Variance estimate of X
## c xp      Plotting Array, Returned
## c         xp(.,1) mean - 1 standard deviation
## c         xp(.,2) smoothed x
## c         xp(.,3) mean + 1 standard deviation
## c nehl    Number each half for smoothing limits
## -- for smoothing variance
## c nehc    Number Each Half for Center
## -- for smoothing data
## c slo     Low end symmetry - "even","odd", "extend"
## c shi   High end symmetry
.lftr3p <- function(x, xv, nhi, nehl, nehc, slo, shi) {
    ip <- 2

    xp <- matrix(NA, nhi, 3)

    xp[,3] <- .llftr7(xv, nhi, "hi", "even", "extend", nehl, ip)
    xp[,2] <- .llftr7(x, nhi, "-", slo, shi, nehc, ip)

    xp[,1] <- xp[,2] - sqrt(xp[,3])
    xp[,3] <- xp[,2] + sqrt(xp[,3])

    return(xp)
}


##modified djt jk
##Programmers note. currently uses jkcoh7 which containes segment averaging
##Segment averaging can be looked into and implemented
##Segment averaging can likely be vectorized or should it be
##implemented in C?
## Dave sets ip=2 in the file ts/lftr3p.f
.llftr7 <- function(x,nhi,lohi,slo,shi,neh,ip) {

    nlo <- 1
    y <- array(NA, nhi)
    z <- array(NA, nhi+2*neh)

    zNlo <- nlo + neh
    zNhi <- nhi + neh
    zHi <- zNhi + neh

    ##  Generate Weights
    fw <- as.double(neh + 1)
    wt <- ((1-((-neh:neh)/fw))*(1+((-neh:neh)/fw)))**ip
    cwt <- sum(wt)
    wt <-  wt/cwt
    
    ## Move data to working (z) array, extend ends,default
    innerSeq <- zNlo:zNhi
    z[innerSeq] <- x
    lowerSeq <- 1:neh
    z[lowerSeq] <- x[1]
    upperSeq <- (zNhi+1):zHi
    z[upperSeq] <- x[nhi]
    
    ## Low End
    
    if(tolower(slo) == "even") {
        z[rev(lowerSeq)] <-  z[(zNlo+1):(zNlo+neh)]
    }

    
    if(tolower(slo) == "odd") {
        ## not tested...
        z[rev(lowerSeq)] <-  2.0*z[zNlo] - z[(zNlo+1):(zNlo+neh)]
    }

    ## High End
    ## Programmer's note: The numbers map to numbers on do loops
    ## in Thomson's fortran code.
    if(tolower(shi) == "even") {
        for( j in  1:neh) { ## 850
            z[zNhi+j] <- z[zNhi-j]
        } ## 850 
    }

    if(tolower(shi) == "odd") {
        for( j  in 1:neh) { ## 860
            z[zNhi + j] =  2.*z[zNhi] - z[zNhi-j]
        } ## 860
    }
    
    ## High Limit, supress local Minima
    if(tolower(lohi) ==  "hi") {
        for( n in  (nlo+1):(nhi-1) ) { ## 1400
            if( (x[n] < x[n-1]) &&  (x[n]  < x[n+1]) ) {
                zNoff <- n +neh
                z[zNoff] <- (x[n-1]+x[n+1])/2.0
            }
        } ## 1400
    }
    ## Low Limit, supress local Maxima
    if(tolower(lohi) == "lo") {
        for( n in  (nlo+1):(nhi-1) ) {
            if( (x[n] > x[n-1] ) &&  (x[n] > x[n+1]) ) {
                zNoff <- n +neh
                z[n] = (x[n-1]+x[n+1])/2.0
            }
        } ## 1500
    }

    zOffSetSeq <- 1:(2*neh+1)
    for (n in nlo:nhi) { ## 2000 
        y[n] <- sum(wt*z[zOffSetSeq])
        zOffSetSeq <- zOffSetSeq +1
    } ## 2000
    return(y)
    ##checks out on first test
    
}


### functions to convert the coherence to the transformed coherence
## normalizing constant 2k-2
.trnrm <- function(k) sqrt(2*k-2)

## These formulae are based on:
## "Jackknifed error estimates for spectra, coherences,
## and transfer functions"
## by Thomson, DJ and Chave, AD
## Advances in Spectrum Estimation
##

## coherence to quantiles of the CDF
.C2toF <-  function(xx, trnrm_) {
    return( trnrm_*log((1.0+sqrt(xx))/(1.0-sqrt(xx)))/2.0 )
}

##quantiles to MSC
.FtoMSC <- function(ff, trnrm_) tanh(ff/trnrm_)**2


## odinlibs nplot function....
## c
## c       Cumulative probability points in 1-2-5 sequence
## c ndata Number data points; output approximately from 1/ndata
## c       to 1 - 1/ndata
## c nmax  Maximum number of Outputs = Dimension of out,cout,Qnorm
## c       nmax approx > 6* log10(ndata)
## c
## c out   Cumulative distribution
## c cout  Character*8 version of CDF
## c Qnorm Quantiles of Standard Normal at CDF
## c nout  Number of output points
## c
## using jkcoh defaults...
.paxpt7 <- function(ndata=2000, nmax=40) {
    ndec <-  round(log10(max(11,ndata)));
    nout = 6*ndec -1;
    if(nout > nmax) {
        return;
    }
    n <- 0;
    out <- array(NA, nout)
    for(m in seq(-ndec, -1, 1)) {
        for(k in c(1,2,5)) {
            n <- n +1;
            v <- as.double(k*10**m);
            out[n] <- v;
            out[nout +1 -n] <- as.double(1.0 - v);
        }
    }
    return(list(out=out, Qnorm=qnorm(out),nout=nout));
}

.cdfToMSqCoh <- function(cdf, k) {
    fnavm <- as.double(k-1);
    return(1.0  - (1.0 - cdf)**(1.0/fnavm));
}

## used in coherence plot--added May, 2013.
.mscToCDFquantiles <- function(msc, k) {
    1 - (1-msc)^(k-1)
}
## end utilities added mainly for plot.mtm.coh
