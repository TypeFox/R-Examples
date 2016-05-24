#
# An implementation of the OPI that simulates responses using 
# Henson et al (2000) variability and also returns response times
# using data from McKednrick et al 2014.
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: August 2013
#
# Modified Tue  8 Jul 2014: added type="X" to opiInitialise and opiPresent
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

simH_RT.opiClose         <- function() { return(NULL) }
simH_RT.opiQueryDevice   <- function() { return (list(type="SimHensonRT")) }

if (!exists(".SimHRTEnv"))
    .SimHRTEnv <- new.env(size=5)

################################################################################
# Input
#   type N|G|C for the three Henson params
#   type X to specify your own A and B values (eg different dB scale)
#   cap  dB value for capping stdev form Henson formula
#   display Dimensions of plot area (-x,+x,-y,+y) to display stim. No display if NULL
#   rtData data.frame with colnames == "Rt", "Dist", "Person"
#
# Side effects if successful:
#   Set .SimHRTEnv$type   to type
#   Set .SimHRTEnv$cap    to cap
#   Set .SimHRTEnv$A      to A
#   Set .SimHRTEnv$B      to B
#   Set .SimHRTEnv$rtData to 3 col data frame to rtData
#
# Return NULL if successful, string error message otherwise  
################################################################################
simH_RT.opiInitialize <- function(type="C", cap=6, A=NA, B=NA, display=NULL, maxStim=10000/pi, rtData, rtFP=1:1600) {
    if (!is.element(type,c("N","G","C", "X"))) {
        msg <- paste("Bad 'type' specified for SimHensonRT in opiInitialize():",type)
        warning(msg)
        return(msg)
    }

    if (cap < 0)
        warning("cap is negative in call to opiInitialize (SimHensonRT)")
    .SimHRTEnv$type <- type
    .SimHRTEnv$cap  <-  cap
    .SimHRTEnv$A    <-  A
    .SimHRTEnv$B    <-  B
    .SimHRTEnv$maxStim <- maxStim

    if (type == "X" && (is.na(A) || is.na(B)))
        warning("opiInitialize (SimHenson): you have chosen type X, but one/both A and B are NA")

    if(simDisplay.setupDisplay(display))
        warning("opiInitialize (SimHensonRT): display parameter may not contain 4 numbers.")

    #if (rtType == "sigma") {
    #    load(paste(.Library,"/OPI/data/RtSigmaUnits.RData",sep=""), envir=.SimHRTEnv)
    #    assign("rtData", .SimHRTEnv$RtSigmaUnits, envir=.SimHRTEnv)
    #} else if (rtType == "db") {
    #    load(paste(.Library,"/OPI/data/RtDbUnits.RData",sep=""), envir=.SimHRTEnv)
    #    assign("rtData", .SimHRTEnv$RtDbUnits, envir=.SimHRTEnv)
    #} else {
    #    msg <- paste("opiInitialize (SimHensonRT): unknown response time data type",rtType)
    #    warning(msg)
    #    return(msg)
    #}

    if (nrow(rtData) < 100) 
        warning("opiInitialize (SimHensonRT): Less than 100 rows in rtData; that's wierd")
    if (ncol(rtData) != 3 || !all(colnames(rtData) == c("Rt", "Dist", "Person"))) {
        msg <- "opiInitialize (SimHensonRT): rtData must have 3 columns: Rt, Dist, Person. See data(RtSigmaUnits) for example."
        warning(msg)
        return(msg)
    }
    assign("rtData", rtData, envir=.SimHRTEnv)

    #print(.SimHRTEnv$rtData[1:10,])

    if (length(rtFP) < 1) {
        msg <- "opiInitialize (SimHensonRT): rtFP must have at least 1 element"
        warning(msg)
        return(msg)
    }

    assign("rtFP", rtFP, envir=.SimHRTEnv)

    return(NULL)
}

################################################################################
# Set background of plot area to col
# Return:
#   NULL - succsess
#   -1   - opiInitialize not called
################################################################################
simH_RT.opiSetBackground <- function(col, gridCol) { 
    return (simDisplay.setBackground(col, gridCol))
}

################################################################################
#
################################################################################
simH_RT.opiPresent <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30) { 
                            UseMethod("simH_RT.opiPresent") }
setGeneric("simH_RT.opiPresent")

#
# Helper function that allows different coefficients from Table 1 of Henson 2000.
# Note prob seeing <0 is always false positive rate (but false neg still poss)
# Response time for false positive is uniform sample from .SimHRTEnv$rtFP
#
simH_RT.present <- function(db, cap=6, fpr=0.03, fnr=0.01, tt=30, dist, A, B) {

    falsePosRt <- function() {
        if(length(.SimHRTEnv$rtFP) < 2) 
            return(.SimHRTEnv$rtFP)
        else
            return(sample(.SimHRTEnv$rtFP,1))
    }

    if (tt < 0)         # force false pos if t < 0
        fpr <- 1.00  

    if (runif(1) < 0.5) {
            # test fp 
        if (runif(1) < 2*fpr) {
            return(list(err=NULL, seen=TRUE, time=falsePosRt()))  # false P
        }
    } else {
            # test fn 
        if (runif(1) < fnr) {
            return(list(err=NULL, seen=FALSE, time=0))                         # false N
        }
    }

        # if get to here then need to check Gaussian
        # and if seen=TRUE need to get a time from .SimHRTEnv$rtData
        # assume pxVar is sigma for RT is in sigma units

    pxVar <- min(cap, exp(A*tt + B)) # variability of patient, henson formula 
    if ( runif(1) < 1 - pnorm(db, mean=tt, sd=pxVar)) {

        o <- head(order(abs(.SimHRTEnv$rtData$Dist - dist)), 100)

        return(list(err=NULL, seen=TRUE, time=sample(.SimHRTEnv$rtData[o,"Rt"], 1)))
    } else {
        return(list(err=NULL, seen=FALSE, time=0))
    }
}#simH_RT.present()

#
# stim is list of type opiStaticStimulus
#
simH_RT.opiPresent.opiStaticStimulus <- function(stim, nextStim=NULL, fpr=0.03, fnr=0.01, tt=30, dist=stim$level - tt) {
    if (!exists("type", envir=.SimHRTEnv)) {
        return ( list(
            err = "opiInitialize(type,cap) was not called before opiPresent()",
            seen= NA,
            time= NA 
        ))
    }

    if (is.null(stim))
        stop("stim is NULL in call to opiPresent (using SimHensonRT, opiStaticStimulus)")

    if (length(tt) != length(fpr))
        warning("In opiPresent (using SimHensonRT), recycling tt or fpr as lengths differ")
    if (length(tt) != length(fnr))
        warning("In opiPresent (using SimHensonRT), recycling tt or fnr as lengths differ")
    if (length(fpr) != length(fnr))
        warning("In opiPresent (using SimHensonRT), recycling fpr or fnr as lengths differ")

    simDisplay.present(stim$x, stim$y, stim$color, stim$duration, stim$responseWindow)

    if (.SimHRTEnv$type == "N") {
        return(simH_RT.present(cdTodb(stim$level, .SimHRTEnv$maxStim), .SimHRTEnv$cap, fpr, fnr, tt, dist, -0.066, 2.81))
    } else if (.SimHRTEnv$type == "G") {
        return(simH_RT.present(cdTodb(stim$level, .SimHRTEnv$maxStim), .SimHRTEnv$cap, fpr, fnr, tt, dist, -0.098, 3.62))
    } else if (.SimHRTEnv$type == "C") {
        return(simH_RT.present(cdTodb(stim$level, .SimHRTEnv$maxStim), .SimHRTEnv$cap, fpr, fnr, tt, dist, -0.081, 3.27))
    } else if (.SimHRTEnv$type == "X") {
        return(simH_RT.present(cdTodb(stim$level, .SimHRTEnv$maxStim), .SimHRTEnv$cap, fpr, fnr, tt, dist, .SimHRTEnv$A, .SimHRTEnv$B))
    } else {
        return ( list(
            err = "Unknown error in opiPresent() for SimHensonRT",
            seen= NA,
            time= NA 
        ))
    }
}

########################################## TO DO !
simH_RT.opiPresent.opiTemporalStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simH_RT temporal persenter yet")
}

simH_RT.opiPresent.opiKineticStimulus <- function(stim, nextStim=NULL, ...) {
    stop("ERROR: haven't written simH_RT kinetic persenter yet")
}

