##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
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
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com


################################################################
##
##  .sphsed
##
##  Phase wrapping routine; takes phases and tracks violation
##  of +/-360 degree boundary, and wraps aliases. For use
##  by demod.dpss().
##
################################################################
.sphsed <-  function(ph,nfreq=length(ph)) {

    q <- 0.0
    pinc <- 0.0

    for(n in 1:nfreq) {
        t1 <- ph[n]
        d <- q - t1
        q <- t1
        if(abs(d) > 180.0) {
            pinc <-  pinc + sign(d)*360.0
        }
        ph[n] <- t1+pinc
    }
    return(ph)
}

################################################################
##
##  demod.dpss
##
##  Complex demodulation routine. Takes a series x, and 
##  demodulates the series around center frequency centreFreq,
##  using parameters NW, blockLen, and stepSize.
##
################################################################
demod.dpss <- function(x, 
                       centreFreq, 
                       nw, 
                       blockLen, 
                       stepSize=1, 
                       wrapphase=TRUE,
                       ...) {

    stopifnot(stepSize == 1)  ## not implemented

    nwTmp <- match.call(expand.dots = )$NW
    
    if(!is.null(nwTmp)) {
        warning("NW has been depreciated. Please use nw instead.")
        nw <- nwTmp
    }
    
    ndata <- length(x)

    deltaT <- deltat(x)
    dw <- dpss(blockLen, 1, nw)$v
    U0 <- sum(dw)

    ampScale <- 2.0/U0
    omegaDeltaT <- 2*pi*centreFreq*deltaT
    jSeq <- (1:blockLen) -1

    complexVal <- exp(-1i*omegaDeltaT*jSeq)
    complexVal <- complexVal*dw*ampScale

    nResultVals <- ndata - blockLen +1

    complexDemod <- complex(nResultVals)

    for(i in 1:nResultVals) {
        iSeq <- i:(i+blockLen-1)
        complexDemod[i] <- crossprod(x[iSeq], complexVal)
    }

    phase <- Arg(complexDemod)*180/pi
    if(wrapphase) {
        phase <- .sphsed(phase)
    }
    phase <- phase - 360*deltaT*centreFreq * (1:nResultVals) 

    list(amplitude=Mod(complexDemod), phase=phase, complexDemod=complexDemod)
}
