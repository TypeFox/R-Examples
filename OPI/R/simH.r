#
# An implementation of the OPI that simulates responses using 
# Henson et al (2000) variability.
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: June 2012
#
# Modified  8 Jul 2014: added type="X" to opiInitialise and opiPresent
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

simH.opiClose         <- function() { return(NULL) }
simH.opiQueryDevice   <- function() { return (list(type="SimHenson")) }

if (!exists(".SimHEnv"))
    .SimHEnv <- new.env(size=5)

################################################################################
# Input
#   type N|G|C for the three Henson params
#   type X to specify your own A and B values (eg different dB scale)
#   cap  dB value for capping stdev form Henson formula
#   display Dimensions of plot area (-x,+x,-y,+y) to display stim. No display if NULL
#   maxStim Maximum stimulus value in cd/m^2 used for db <-> cd/m^2 conversions
#
# Return NULL if succesful, string error message otherwise  
################################################################################
simH.opiInitialize <- function(type="C", A=NA, B=NA, cap=6, display=NULL, maxStim=10000/pi) {
    if (!is.element(type,c("N","G","C","X"))) {
        msg <- paste("Bad 'type' specified for SimHenson in opiInitialize():",type)
        warning(msg)
        return(msg)
    }

    if (cap < 0)
        warning("cap is negative in call to opiInitialize (simHenson)")
    .SimHEnv$type <- type
    .SimHEnv$cap  <-  cap
    .SimHEnv$A    <-  A
    .SimHEnv$B    <-  B
    .SimHEnv$maxStim <- maxStim

    if (type == "X" && (is.na(A) || is.na(B)))
        warning("opiInitialize (SimHenson): you have chosen type X, but one/both A and B are NA")

    if(simDisplay.setupDisplay(display))
        warning("opiInitialize (SimHenson): display parameter may not contain 4 numbers.")

    return(NULL)
}

################################################################################
# Set background of plot area to col
# Return:
#   NULL - succsess
#   -1   - opiInitialize not called
################################################################################
simH.opiSetBackground <- function(col, gridCol) { 
    return (simDisplay.setBackground(col, gridCol))
}

################################################################################
#
################################################################################
simH.opiPresent <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30) { 
                            UseMethod("simH.opiPresent") }
setGeneric("simH.opiPresent")

#
# Helper function that allows different coefficients from Table 1 of Henson 2000.
# Note prob seeing <0 is always false positive rate 
#
simH.present <- function(db, cap=6, fpr=0.03, fnr=0.01, tt=30, A, B) {

    if (tt >= 0) {
        pxVar <- min(cap, exp(A*tt + B)) # variability of patient, henson formula 

        prSeeing <- fpr + (1-fpr-fnr)*(1-pnorm(db, mean=tt, sd=pxVar))    
    } else {
        prSeeing <- fpr
    }

    return ( list(
        err = NULL,
        seen= runif(1) < prSeeing,
        time= 0
    ))
}#

#
# stim is list of type opiStaticStimulus
#
simH.opiPresent.opiStaticStimulus <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30) {
    if (!exists("type", envir=.SimHEnv)) {
        return ( list(
            err = "opiInitialize(type,cap) was not called before opiPresent()",
            seen= NA,
            time= NA 
        ))
    }

    if (is.null(stim))
        stop("stim is NULL in call to opiPresent (using simHenson, opiStaticStimulus)")

    if (length(tt) != length(fpr))
        warning("In opiPresent (using simHenson), recycling tt or fpr as lengths differ")
    if (length(tt) != length(fnr))
        warning("In opiPresent (using simHenson), recycling tt or fnr as lengths differ")
    if (length(fpr) != length(fnr))
        warning("In opiPresent (using simHenson), recycling fpr or fnr as lengths differ")

    simDisplay.present(stim$x, stim$y, stim$color, stim$duration, stim$responseWindow)

    if (.SimHEnv$type == "N") {
        return(simH.present(cdTodb(stim$level, .SimHEnv$maxStim), .SimHEnv$cap, fpr, fnr, tt, -0.066, 2.81))
    } else if (.SimHEnv$type == "G") {
        return(simH.present(cdTodb(stim$level, .SimHEnv$maxStim), .SimHEnv$cap, fpr, fnr, tt, -0.098, 3.62))
    } else if (.SimHEnv$type == "C") {
        return(simH.present(cdTodb(stim$level, .SimHEnv$maxStim), .SimHEnv$cap, fpr, fnr, tt, -0.081, 3.27))
    } else if (.SimHEnv$type == "X") {
        return(simH.present(cdTodb(stim$level, .SimHEnv$maxStim), .SimHEnv$cap, fpr, fnr, tt, .SimHEnv$A, .SimHEnv$B))
    } else {
        return ( list(
            err = "Unknown error in opiPresent() for SimHenson",
            seen= NA,
            time= NA 
        ))
    }
}

########################################## TO DO !
simH.opiPresent.opiTemporalStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simH temporal persenter yet")
}

simH.opiPresent.opiKineticStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simH kinetic persenter yet")
}

