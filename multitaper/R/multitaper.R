##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Karim Rahim and Wesley S. Burr.
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
##     

##############################################################
##
##  spec.mtm
##
##  Wrapper routine for .spec.mtm.dpss and .spec.mtm.sine.
##
##############################################################
spec.mtm <- function(timeSeries,
                     nw=4.0,
                     k=7,
                     nFFT="default", 
                     taper=c("dpss"),
                     centre=c("Slepian"),
                     dpssIN=NULL,
                     returnZeroFreq=TRUE,
                     Ftest=FALSE,
                     jackknife=FALSE,
                     jkCIProb=.95,
                     adaptiveWeighting=TRUE,
                     maxAdaptiveIterations=100,
                     plot=TRUE,
                     na.action=na.fail,
                     returnInternals=FALSE,
                     sineAdaptive=FALSE,
                     sineSmoothFact=0.2,
                     dtUnits=c("default"),
                     deltat=NULL,
                     ...) {
    
    series <- deparse(substitute(timeSeries))
    taper <- match.arg(taper,c("dpss","sine"))
    centre <- match.arg(centre,c("Slepian","arithMean","trimMean","none"))
    dtUnits <- match.arg(dtUnits,c("second","hour","day","month","year","default"))

    ## deal with depreciated parameter dT is changed to deltat
    ## we strip dT before plotting in plotsHelper.R
    ## to prevent it getting passed to plot
 
    deltaT <- NULL
    if(!missing(deltat)) {
        deltaT <- deltat
    }
   
    dT <- match.call(expand.dots = )$dT
    
    if(missing(deltat) && !is.null(dT)) {
        warning("dT has been depreciated. Use either deltat or input a time series object.")
        deltaT <- dT
    }
    
    if( (taper=="sine") && is.complex(timeSeries)) {
        stop("Sine tapering not implemented for complex time series.") 
    }
    if( (taper=="sine") && jackknife) { 
        warning("Cannot jackknife over sine tapers.")
        jackknife <- FALSE
    }
    if( (taper=="sine") && Ftest) { 
        warning("Cannot compute Ftest over sine tapers.")
        Ftest <- FALSE
    } 
    if( (taper=="sine") && !returnZeroFreq) {
        returnZeroFreq = TRUE 
        warning("returnZeroFreq must be TRUE for sine taper option.")
    }
    if( (taper=="sine") && sineSmoothFact > 0.5) {
        warning("Smoothing Factor > 0.5 is very high!")
    }
    
    dtTmp <- NULL
    ## warning for deltaT missing: makes all frequency plots incorrect
    if(!is.ts(timeSeries) && is.null(deltaT)) {
        warning("Time series is not a ts object and deltat is not set. Frequency array and axes may be incorrect.")
    }
    if(!is.ts(timeSeries)) {
        if(!is.complex(timeSeries)) {
            timeSeries <- as.double(as.ts(timeSeries)) 
        }
    } else {
        ## Order matters here, because as.double breaks the ts() class
        dtTmp <- deltat(timeSeries)
        if(!is.complex(timeSeries)) {
            timeSeries <- as.double(timeSeries)
        }
    }
    
    ## in responese to delta T bug July 2, 2013
    ## modified to remove dT
    
    if(is.null(deltaT)) {
        if(!is.null(dtTmp)) {
            deltaT <- dtTmp
        } else{
            deltaT <- 1.0
        }           
    }
    
    n <- length(timeSeries)

    if(taper=="dpss") {
        stopifnot(nw >= 0.5, k >= 1, n > 8)
        ## replace stop if not with warning.
        ## the following was also in stopif not:
        ## nw <= 500, k <= 1.5+2*nw)
        if( nw > 500) {
            warning("nw > 500")
        }
        if( k > 1.5 * 2*nw ) {
            warning("k > 1.5+2*nw")
        }
        if (nw/n > 0.5) { 
            warning("half-bandwidth parameter (w) is greater than 1/2")
        }
        if(k==1) {
            Ftest=FALSE
            jackknife=FALSE
        }
    } else {
        stopifnot(k <= n, k >= 1, n > 8)
    }

    na.action(timeSeries)
    if(!is.complex(timeSeries)) {
        sigma2 <- var(timeSeries) * (n-1)/n
    } else {
        sigma2 <- var(Re(timeSeries)) * (n-1)/n + var(Im(timeSeries)) * (n-1)/n 
    }

    if(nFFT == "default") {
        nFFT <- 2* 2^ceiling(log2(n))
    } else {
        stopifnot(is.numeric(nFFT))
    }
    stopifnot(nFFT >= n)
    
    ## convert time-series to zero-mean by one of three methods, if set; default is Slepian
    if(centre=="Slepian") {
        if(taper=="dpss") {
            timeSeries <- centre(timeSeries, nw=nw, k=k, deltaT=deltaT)    
        } else {  # edge case: sine taper, set initial k, but too high for default nw=4.0
            timeSeries <- centre(timeSeries, nw=5.0, k=8, deltaT=deltaT)
        }
    } else if(centre=="arithMean") {
        timeSeries <- centre(timeSeries, trim=0) 
    } else if(centre=="trimMean") {
        timeSeries <- centre(timeSeries, trim=0.10)
    }
   
    if(taper=="dpss") { 
        mtm.obj <- .spec.mtm.dpss(timeSeries=timeSeries,
                                  nw=nw, k=k, nFFT=nFFT, 
                                  dpssIN=dpssIN, returnZeroFreq=returnZeroFreq, 
                                  Ftest=Ftest, jackknife=jackknife, jkCIProb=jkCIProb, 
                                  adaptiveWeighting = adaptiveWeighting, 
                                  maxAdaptiveIterations=maxAdaptiveIterations, 
                                  returnInternals=returnInternals, 
                                  n=n, deltaT=deltaT, sigma2=sigma2, series=series,
                                  dtUnits=dtUnits, ...) 
    } else if(taper=="sine") {
        mtm.obj <- .spec.mtm.sine(timeSeries=timeSeries, k=k, sineAdaptive=sineAdaptive,
                                  nFFT=nFFT, dpssIN=dpssIN, returnZeroFreq=returnZeroFreq,
                                  returnInternals=FALSE, n=n, deltaT=deltaT, sigma2=sigma2,
                                  series=series,maxAdaptiveIterations=maxAdaptiveIterations,
                                  smoothFact=sineSmoothFact, dtUnits=dtUnits, ...)
    }

    if(plot) {
        plot.mtm(mtm.obj, jackknife=jackknife, ...)
        return(invisible(mtm.obj))
    } else {
        return(mtm.obj)
    }
}

##############################################################
##
##  .spec.mtm.dpss
##
##  Computes multitaper spectrum using Slepian tapers
##  References: 
##    Percival and Walden "Spectral Analysis
##    for Physical Applications" 1993 and associated LISP code
##
##    Thomson, D.J. Spectrum Estimation and Harmonic Analysis,
##    Proceedings of the IEEE, 1982 and associated Fortran code
## 
##############################################################
.spec.mtm.dpss <- function(timeSeries,
                     nw,
                     k,
                     nFFT,
                     dpssIN,
                     returnZeroFreq,
                     Ftest,
                     jackknife,
                     jkCIProb,
                     adaptiveWeighting, 
                     maxAdaptiveIterations,
                     returnInternals,
                     n,
                     deltaT,
                     sigma2,
                     series,
                     dtUnits,
                     ...) {

    # Complex check case
    if(is.complex(timeSeries)) {
      if(!returnZeroFreq) {
        returnZeroFreq <- 1 
        warning("Cannot set returnZeroFreq to 0 for complex time series.")
      } 
    }

    dw <- NULL
    ev <- NULL
    receivedDW <- TRUE

    if(!.is.dpss(dpssIN)) {
      receivedDW <- FALSE
      dpssIN <- dpss(n, k, nw=nw, returnEigenvalues=TRUE)
      dw <- dpssIN$v*sqrt(deltaT)
      ev <- dpssIN$eigen 
    }
    else {
      dw <- .dpssV(dpssIN)
      ev <- .dpssEigen(dpssIN)
      if(all(is.null(ev))) {
        ev <- dpssToEigenvalues(dw, nw) }
        dw <- dw*sqrt(deltaT) 
    }

    nFreqs <- nFFT %/% 2 + as.numeric(returnZeroFreq)
    offSet <- if(returnZeroFreq) 0 else 1 

    # Note that the frequency axis is set by default to unit-less
    # scaling as 0 through 0.5 cycles/period. The user parameter
    # dtUnits modifies this scaling in the plot.mtm function.
    scaleFreq <- 1 / as.double(nFFT * deltaT)
    
    swz <- NULL ## Percival and Walden H0
    ssqswz <- NULL
    swz <- apply(dw, 2, sum)
    if(k >= 2) {
      swz[seq(2,k,2)] <- 0
    }
    ssqswz <- as.numeric(t(swz)%*%swz)

    taperedData <- dw * timeSeries
    
    nPadLen <- nFFT - n
    if(!is.complex(timeSeries)) {
      paddedTaperedData <- rbind(taperedData, matrix(0, nPadLen, k))
    } else {
      paddedTaperedData <- rbind(taperedData, matrix(complex(0,0), nPadLen, k)) 
    }
    cft <- mvfft(paddedTaperedData)
    if(!is.complex(timeSeries)) {
      cft <- cft[(1+offSet):(nFreqs+offSet),]
    } else {
      cft <- rbind(cft[(nFreqs+offSet+1):nFFT,],cft[(1+offSet):(nFreqs+offSet),])
    }
    sa <- abs(cft)^2
   
    if(!is.complex(timeSeries)) {
      resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*scaleFreq 
    } else {
      resultFreqs <- (-(nFreqs-1):(nFreqs-2))*scaleFreq
    }

    adaptive <-  NULL
    jk <- NULL
    PWdofs <- NULL
    if(!jackknife) {
        if(!is.complex(timeSeries)) {
          adaptive <- .mw2wta(sa, nFreqs, k, sigma2, deltaT, ev)
        } else {
          adaptive <- .mw2wta(sa, nFFT, k, sigma2, deltaT, ev) 
        }
    } else {
        stopifnot(jkCIProb < 1, jkCIProb > .5)
        if(!is.complex(timeSeries) & adaptiveWeighting) {
          adaptive <- .mw2jkw(sa, nFreqs, k, sigma2, deltaT, ev)
        } else {
          adaptive <- .mw2jkw(sa, nFFT, k, sigma2, deltaT, ev)
        }
        scl <- exp(qt(jkCIProb,adaptive$dofs)*
                   sqrt(adaptive$varjk))
        upperCI <- adaptive$s*scl
        lowerCI <- adaptive$s/scl
        minVal = min(lowerCI)
        maxVal = max(upperCI)
        jk <- list(varjk=adaptive$varjk,
                   bcjk=adaptive$bcjk,
                   sjk=adaptive$sjk,
                   upperCI=upperCI,
                   lowerCI=lowerCI,
                   maxVal=maxVal,
                   minVal=minVal)
   } 

   ftestRes <- NULL

   if(Ftest) {
        if(is.null(swz)) {
            swz <- apply(dw, 2, sum)
        }
        ftestRes <- .HF4mp1(cft,
                            swz,
                            k,
                            ssqswz)
    }

    eigenCoef1 <- NULL
    wtCoef1 <- NULL
    
    if(returnInternals) {
        eigenCoef1 <- cft
        if(adaptiveWeighting) {
          wtCoef1 <- sqrt(adaptive$wt)
        } else {
          wtCoef <- rep(1, nFreqs)
        }
    }
    auxiliary <- list(dpss=dpssIN,
                      eigenCoefs=eigenCoef1,
                      eigenCoefWt=wtCoef1,
                      nfreqs=nFreqs,
                      nFFT=nFFT,
                      jk=jk,
                      Ftest=ftestRes$Ftest,
                      cmv=ftestRes$cmv,
                      dofs=adaptive$dofs,
                      nw=nw,
                      k=k,
                      deltaT=deltaT,
                      dtUnits=dtUnits,
                      taper="dpss")

    ##   Thomson, D.J. Spectrum Estimation and Harmonic Analysis,
    ##   Proceedings of the IEEE, 1982.

    ## note that the weights are squared, they are |d_k(f)^2 from equation
    ## (5.4)
    ## These weights correspond to Thomoson's 1982 Fortran code.
    ## dof fix for one taper, only value.
    if(k==1) {
        auxiliary$dofs <- 2
    }
    
    spec.out <- list(origin.n=n,
                     method="Multitaper Spectral Estimate",
                     pad= nFFT - n,
                     spec=adaptive$s,
                     freq=resultFreqs,
                     series=series,
                     adaptive=adaptiveWeighting, 
                     mtm=auxiliary)

    class(spec.out) <- c("mtm", "spec")
    
    if(Ftest) {
        class(spec.out) <- c("mtm", "Ftest", "spec")
    }
    return(spec.out)
}


#########################################################################
##
##  spec.mtm.sine 
##
##  Computes multitaper spectrum estimate using sine tapers, as in
## 
##  Riedel, Kurt S. and Sidorenko, Alexander, Minimum Bias Multiple 
##    Taper Spectral Estimation. IEEE Transactions on Signal Processing,
##    Vol. 43, No. 1, January 1995.
##
##  Algorithm implementation based on previous work by:
##    German Prieto, Universidad de los Andes
##       via \texttt{mtsepc}, a F90 package that can be found at
##       http://wwwprof.uniandes.edu.co/~gprieto/software/mwlib.html
##
##    and
##
##    Robert L. Parker, Scripps Institution of Oceanography
##      via \texttt{psd.f}, a F77 program that can be found at
##      http://igppweb.ucsd.edu/~parker/Software/Source/psd.f
## 
#########################################################################

.spec.mtm.sine <- function(timeSeries,
                          nFFT,
                          k,
                          sineAdaptive,
                          dpssIN, 
                          returnZeroFreq=TRUE,
                          n, 
                          deltaT, 
                          dtUnits,
                          sigma2,
                          series=series,
                          maxAdaptiveIterations,
                          smoothFact,
                          ...) {

    dw <- NULL
    receivedDW <- TRUE
    if(!.is.dpss(dpssIN)) {
      receivedDW <- FALSE
      dpssIN <- sineTaper(n, k)
      dw <- dpssIN$v
    }
    else {
      dw <- .dpssV(dpss)
    }

    # returnZeroFreq forced to TRUE, offset = 0
    # NOTE: sine tapers produce nFFT/4 unique results; need to scale nFFT and nFreqs accordingly
    nFFT <- nFFT*2
    nFreqs <- nFFT %/% 4 + as.numeric(returnZeroFreq)
    offSet <- if(returnZeroFreq) 0 else 1 
    scaleFreq <- 1 / as.double(nFFT/2 * deltaT)
    resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*scaleFreq 
    nPadLen <- nFFT - n
    df <- 1/as.double(nFFT*deltaT)

    # compute a single FFT; since we are using sine tapers, this is all we need
    ones <- matrix(1,n,1)
    paddedData<- rbind(timeSeries*ones, matrix(0, nPadLen, 1))
    cft <- mvfft(paddedData)
 
    # constant number of tapers, or adaptive?
    spec <- as.double(matrix(0,1,nFreqs))

    if(!sineAdaptive) { # constant k tapers
       spec <- (.qsF(nFreqs=nFreqs,nFFT=nFFT,k=k,cft=cft,useAdapt=FALSE,kadapt=c(1)))$spec
       dofs <- NULL
    } else { # adaptively weighted tapers

      initTaper <- ceiling(3.0 + sqrt(smoothFact*n)/5.0);

      # pilot estimate of S
      spec0 <- (.qsF(nFreqs=nFreqs,nFFT=nFFT,k=k,cft=cft,useAdapt=FALSE,kadapt=c(1)))$spec

      out <- .adaptSine (ntimes=maxAdaptiveIterations, 
                          k=initTaper, 
                          nFreqs=nFreqs, 
                          sx=spec0, 
                          nFFT=nFFT, 
                          cft=cft,
                          df=df,
                          fact=smoothFact) 
      spec <- out$spec;
      dofs <- out$kadapt;
    } # end of adaptive logic

    # normalize spectrum
    const <- var(timeSeries)/sum(spec)/df
    specFinal <- const*spec

    ## set up return object

    if(sineAdaptive) { 
      method = "Sine-Taper Multitaper Spectrum (adaptive)"
    } else {
      method = paste("Sine-Taper Multitaper Spectrum (k=",k,")",sep="")
    }
    

    auxiliary <- list(dpss=dpssIN,
                      eigenCoefs=NULL,
                      eigenCoefWt=NULL,
                      nfreqs=nFreqs,
                      nFFT=nFFT,
                      jk=NULL,
                      Ftest=NULL,
                      cmv=NULL,
                      dofs=dofs,
                      nw=NULL,
                      k=k,
                      deltaT=deltaT,
                      dtUnits=dtUnits,
                      taper="sine")

    spec.out <- list(origin.n=n,
                     method=method,
                     pad= nFFT - n,
                     spec=specFinal,
                     spec = NULL,
                     freq=resultFreqs,
                     series=series,
                     mtm=auxiliary)

    class(spec.out) <- c("mtm", "spec")
    return(spec.out)
}

#########################################################################
##
## centre
##
## Takes a time series and converts to zero-mean using one of three 
## methods: Slepian projection, arithmetic mean, or trimmed mean.
## 
#########################################################################

centre <- function(x, nw=NULL, k=NULL, deltaT=NULL, trim=0) {
    na.fail(x)
    res <- NULL
    if(is.null(nw) && is.null(k) ) {
        res <- x - mean(x, trim=trim)
    } else {
        if(trim != 0) {
            warning(paste("Ignoring trim =", trim))
        }
        stopifnot(nw >= 0.5, k >= 1, nw <= 500, k <= 1.5+2*nw)
        if (nw/length(x) > 0.5) { 
            stop("half-bandwidth parameter (w) is greater than 1/2")
        }
        if(is.null(deltaT)) {
            if(is.ts(x)) {
                deltaT <- deltat(ts)
            } else {
                warning("deltaT not specified; using deltaT=1.")
                deltaT <- 1.0
            }
        }
        n <- length(x)
        dpssRes <- dpss(n, k=k, nw=nw,
                        returnEigenvalues=TRUE)
        dw <- dpssRes$v*sqrt(deltaT)
        ev <- dpssRes$eigen
        swz <- apply(dw, 2, sum)
        ## zero swz where theoretically zero; odd tapers
        if(k >=2) {
          swz[seq(2,k,2)] <- 0.0
        }
        ssqswz <- sum(swz^2)
        if(!is.complex(x)) {
          res <- .mweave(x, dw, swz,
                         n, k, ssqswz, deltaT)
          res <- x - res$cntr
        } else {
          res.r <- .mweave(Re(x), dw, swz,
                           n, k, ssqswz, deltaT)
          res.i <- .mweave(Im(x), dw, swz,
                           n, k, ssqswz, deltaT)
          res <- x - complex(real=res.r$cntr, imaginary=res.i$cntr)
        }
    }
    return(res)
}


#########################################################################
##
## jackknife coherence and helper smoother and plotting functions
## 
## Example: 
## jkRes <- jkcoh1(r1$auxiliary$cft, r2$auxiliary$cft,
##                 4,2048,4,4096,395)
## pGreater <-  percentjkMSCGreaterThan(jkRes$msc, 4)
## plotJkcoh1(r1$freqs, jkRes$TRmsc, jkRes$NTvar, 4, pGreater)
##
#########################################################################

mtm.coh <- function(mtm1, mtm2, fr=NULL, tau=0, phcorr = TRUE, 
                    plot=TRUE,...) {

    ## note Dave saves the cft
    ## in ./odinlibs-1.1/src/mw/mw2pakt as weighted
    ## 1000 blkcft(n,k,curblk,curset) =
    ##  cft(n*ndecfr,k)*sqrt(wt(n*ndecfr,k))

    ## we require auxiliary data
    if(is.null(mtm1$mtm$eigenCoefs) || is.null(mtm2$mtm$eigenCoefs)) {
        stop("Both mtm objects must have been computed with returnInternals=TRUE.")
    }

    if(mtm1$mtm$k != mtm1$mtm$k) {
        stop("Both mtm objects must have the same value for k.")
    }
    ##k <- mtm1$auxiliary$

    if(mtm1$mtm$nfreqs != mtm1$mtm$nfreqs) {
        stop("Both mtm objects must have the same value for nFFT.")
    }

    nord <- mtm1$mtm$k
    nfreqs <- mtm1$mtm$nfreqs
    cft1 <- mtm1$mtm$eigenCoefs
    cft2 <- mtm2$mtm$eigenCoefs
    
    fr <-  if(is.null(fr))  array(as.double(0), nfreqs) else fr
    
    blklof <-  if(nfreqs %%2 ==0) 1 else 0
    blkhif <- nfreqs -1 + blklof

    nordP2 <- nord +2
    out <- .Fortran("jkcoh1", cft1=as.complex(cft1),
                    cft2=as.complex(cft2), nord=as.integer(nord),
                    blklof=as.integer(blklof), blkhif=as.integer(blkhif),
                    fr=as.double(fr),  tau=as.double(tau),
                    phcorr=as.integer(phcorr),
                    NTmsc=double(nfreqs), NTvar=double(nfreqs),
                    msc=double(nfreqs), ph=double(nfreqs),
                    phvar=double(nfreqs),
                    s1=double(nordP2), s2=double(nordP2),
                    jkmsc=double(nordP2), TRmsc=double(nordP2),
                    bias=double(nfreqs),
                    cx=complex(nordP2),
                    PACKAGE="multitaper")

    auxiliary <- list(nfreqs=mtm1$mtm$nFreqs,
                      nFFT=mtm1$mtm$nFFT,
                      nw=mtm1$mtm$nw,
                      k=mtm1$mtm$k,
                      deltaT=mtm1$mtm$deltaT,
                      dtUnits=mtm1$mtm$dtUnits,
                      taper=mtm1$mtm$taper
                      )


    coh.out <- list(NTmsc=out$NTmsc, NTvar=out$NTvar,
                    msc=out$msc, nfreqs=mtm1$mtm$nfreqs,
                    freq=mtm1$freq, k=nord,
                    ph=out$ph, phvar=out$phvar, mtm=auxiliary)
    class(coh.out) <- "mtm.coh"
    
    
   if(plot) {
        plot.mtm.coh(coh.out, ...)
        return(invisible(coh.out))
    } else {
        return(coh.out)
    }
}


