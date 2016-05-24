#
# An implementation of the OPI that simulates responses using 
# Cummulative Gaussian variability with supplied standard deviation.
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: June 2012
#
# Modified 20 Jul 2014: added maxStim argument for cdTodB conversion
#
# Copyright 2012 Andrew Turpin
# This program is part of the OPI (http://perimetry.org/OPI).
# OPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

simG.opiClose         <- function() { return(NULL) }
simG.opiQueryDevice   <- function() { return (list(type="SimGaussian")) }

if (!exists(".SimGEnv"))
    .SimGEnv <- new.env(size=2)

################################################################################
# Input
#   sd standard deviation for the Gaussian
#   display Dimensions of plot area to display stim. No display if NULL
#
# Return NULL if succesful, string error message otherwise  
################################################################################
simG.opiInitialize <- function(sd, display=NULL, maxStim=10000/pi) {
    if (!is.numeric(sd) || (sd < 0)) {
        msg <- paste("Invalid standard deviation in opiInitialize for SimGaussian:",sd)
        warning(msg)
        return(msg)
    }

    .SimGEnv$sd <- sd
    .SimGEnv$maxStim <- maxStim

    if (simDisplay.setupDisplay(display))
        warning("opiInitialize (SimGaussian): perhaps display parameter does not contain 4 numbers?")

    return(NULL)
}

################################################################################
# Set background of plot area to col, color of gird lines to gridCol
################################################################################
simG.opiSetBackground <- function(col, gridCol) { 
    simDisplay.setBackground(col, gridCol)
    return(NULL) 
}
################################################################################
simG.opiPresent <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30) { UseMethod("simG.opiPresent") }
setGeneric("simG.opiPresent")

simG.opiPresent.opiStaticStimulus <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30) {
    if (!exists("sd", envir=.SimGEnv)) {
        return ( list(
            err = "opiInitialize(sd) was not called before opiPresent()",
            seen= NA,
            time= NA 
        ))
    }

    if (is.null(stim))
        stop("stim is NULL in call to opiPresent (using simGaussian, opiStaticStimulus)")

    if (as.numeric(.SimGEnv$sd) <= 0)
        warning("sd is <= 0 in SimGaussian call to opiPresent")

    prSeeing <- fpr + (1-fpr-fnr)*(1-pnorm(cdTodb(stim$level, .SimGEnv$maxStim), mean=tt, sd=as.numeric(.SimGEnv$sd)))

    simDisplay.present(stim$x, stim$y, stim$color, stim$duration, stim$responseWindow)

    return ( list(
        err = NULL,
        seen= runif(1) < prSeeing,
        time= 0
    ))
}#

########################################## TO DO
simG.opiPresent.opiTemporalStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simG temporal persenter yet")
}

simG.opiPresent.opiKineticStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simG kinetic persenter yet")
}
