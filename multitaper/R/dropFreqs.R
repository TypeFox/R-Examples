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
##  dropFreqs.*
##
##  Plotting utility functions that allow the user to subselect
##  a frequency range of interest, and 'drop' the extraneous
##  frequencies. Note that these functions are intended to be
##  used only at the end of analysis, as once they have been
##  applied to an object, the result is not suitable for 
##  passing into any further computational routines (such as 
##  mtm.coh).
##
################################################################


# Handler
dropFreqs <- function(spec, minFreq, maxFreq) UseMethod("dropFreqs")

# Fall-through case
dropFreqs.default <- function(spec, minFreq, maxFreq) {
    print("This function is only valid for objects of spec, mtm, or mtm.coh classes")
    spec
}

# Spectrum object
dropFreqs.spec <- function(spec, minFreq, maxFreq) {
    idx <- (findInterval(spec$freq, c(minFreq,maxFreq)) == 1)

    if(sum(idx) <= 1) {
        stop("minFreq and maxFreq must allow for a range of frequencies to be returned")
    }
    spec.out <- spec
    spec.out$freq <- spec$freq[idx]
    spec.out$spec <- spec$spec[idx]

    spec.out
}

# mtm object
dropFreqs.mtm <- function(spec, minFreq, maxFreq) {
    idx <- (findInterval(spec$freq, c(minFreq,maxFreq)) == 1)

    if(sum(idx) <= 1) {
        stop("minFreq and maxFreq must allow for a range of frequencies to be returned")
    }
    spec.out <- spec
    spec.out$freq <- spec$freq[idx]
    spec.out$spec <- spec$spec[idx]

    ##adjust mtm parameters
    if(!is.null(spec.out$mtm)) {
        ## null unnecessary values 
        ## enforces fact that currently function is mainly a
        ## plotting utility
        spec.out$mtm$dpss <- NULL
        spec.out$mtm$eigenCoefs <- NULL
        spec.out$mtm$eigenCoefsWt <- NULL
        
        ## keep values used in plotting
        spec.out$mtm$Ftest <- spec.out$mtm$Ftest[idx]
        spec.out$mtm$dofs <- spec.out$mtm$dofs[idx]

        if(!is.null(spec.out$mtm$jk)) {
            spec.out$mtm$jk$varjk <- NULL
            spec.out$mtm$jk$upperCI <- spec.out$mtm$jk$upperCI[idx]
            spec.out$mtm$jk$maxVal <- max(spec.out$mtm$jk$upperCI)
            spec.out$mtm$jk$bcjk <- NULL
            spec.out$mtm$jk$lowerCI <- spec.out$mtm$jk$lowerCI[idx]
            spec.out$mtm$jk$sjk <- NULL
            spec.out$mtm$jk$minVal <- min(spec.out$mtm$jk$lowerCI)
        }
    }
    spec.out
}

# mtm.coh object
dropFreqs.mtm.coh <- function(spec, minFreq, maxFreq) {
    idx <- (findInterval(spec$freq, c(minFreq,maxFreq)) == 1)

    if(sum(idx) <= 1) {
        stop("minFreq and maxFreq must allow for a range of frequencies to be returned")
    }

    spec.out <- spec
    spec.out$NTmsc <- spec.out$NTmsc[idx]
    spec.out$msc <- spec.out$msc[idx]
    spec.out$NTvar <- spec.out$NTvar[idx]
    spec.out$freq <- spec.out$freq[idx]
    spec.out$ph <- spec.out$ph[idx]
    spec.out$phvar <- spec.out$phvar[idx]
    spec.out$nfreqs <- sum(idx)
    spec.out
}
