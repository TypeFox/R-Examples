#
# Full Threshold (FT) algorithm that maintains state in a list.
# See comments in full_threshold.r.
# Includes
#     FT.start    # initialise list state
#     FT.step     # take state, present stim, update and return state
#     FT.stop     # boolean - true if state is finished
#     FT.final    # return final estimate from state
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: September 2013
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

################################################################################
# Input parameters
#   est           Starting estimate in dB
#   instRange     Dynamic range of the instrument c(min,max) in dB
#   makeStim      A helper function to create the required
#                 OPI data type for passing to opiPresent
#   ...           Parameters for opiPresent
# Returns a list containing
#   npres    Total number of presentations
#   respSeq  Response sequence stored as a list of (seen,dB) pairs
#   first    First staircase estimate in dB
#   final    Final threshold estimate in dB
################################################################################
FT.start <- function(est=25, instRange=c(0,40), makeStim, ...) {
    if (est < instRange[1] || est > instRange[2])
        stop("FT.start: est must be in the range of instRange")

    return(list(name="FT",
        startingEstimate=est,
        currentLevel=est,
        minStimulus=instRange[1],
        maxStimulus=instRange[2],
        makeStim=makeStim,
        lastSeen=NA,
        lastResponse=NA,
        firstStairResult=NA,
        secondStairResult=NA,
        finished=FALSE,                     # flag to say it is finished
        numberOfReversals=0,
        currSeenLimit=0,                    # number of times maxStimulus seen
        currNotSeenLimit=0,                 # number of times minStimulus not seen
        numPresentations=0,                 # number of presentations so far
        stimuli=NULL,                       # vector of stims shown
        responses=NULL,                     # vector of responses (1 seen, 0 not)
        responseTimes=NULL,                 # vector of response times
        opiParams=list(...)                 # the extra params
    ))
}# FT.start()

################################################################################
# Present state$currentLevel and update state after response.
#
# Input parameters
#   State list as returned by FT.start
#   nextStim - suitable for passing to opiPresent, can be NULL
# Returns a list containing
#   state = updated state list
#   resp  = response to the presentation (as returned by opiPresent)
#
# Follows FT algorithm as documented in full_threshold.r.
# Note 
#   1) opiPresent called infinitely until no error
################################################################################
FT.step <- function(state, nextStim=NULL) {
    if (is.null(state$opiParams))
        params <- list(stim=state$makeStim(state$currentLevel, state$numPresentations), nextStim=nextStim)
    else
        params <- c(list(stim=state$makeStim(state$currentLevel, state$numPresentations), nextStim=nextStim), state$opiParams)
    opiResp <- do.call(opiPresent, params)
    while (!is.null(opiResp$err))
        opiResp <- do.call(opiPresent, params)

    state$stimuli          <- c(state$stimuli, state$currentLevel)
    state$responses        <- c(state$responses, opiResp$seen)
    state$responseTimes    <- c(state$responseTimes, opiResp$time)
    state$numPresentations <- state$numPresentations + 1

    if (opiResp$seen) 
        state$lastSeen <- state$currentLevel

        # check for seeing min
    if (state$currentLevel == state$minStimulus && !opiResp$seen)
        state$currNotSeenLimit <- state$currNotSeenLimit + 1

        # check for seeing max
    if (state$currentLevel == state$maxStimulus && opiResp$seen)
        state$currSeenLimit <- state$currSeenLimit + 1

        # check for reversals
    if (state$numPresentations > 1 && opiResp$seen != state$lastResponse)
        state$numberOfReversals <- state$numberOfReversals + 1 

    state$lastResponse <- opiResp$seen

        # check if staircase is finished.
        # If it is, and it is the first, 
        # check if we need a second staircase
    thisStairStops <- state$numberOfReversals >= 2 || state$currNotSeenLimit >= 2 || state$currSeenLimit >= 2
    if (!thisStairStops) {
            # update stimulus for next presentation
        delta <- ifelse(state$numberOfReversals == 0, 4, 2) * ifelse(opiResp$seen, +1, -1)
        state$currentLevel <- min(state$maxStimulus, max(state$minStimulus, state$currentLevel + delta))
    } else {
        if (is.na(state$firstStairResult)) {
            first <- FT.final.details(state)
            state$firstStairResult <- first$final
            if (first$stopReason == "Reversals" && abs(state$firstStairResult - state$startingEstimate) > 4) { 
                    # initiate second staircase
                state$currentLevel      <- first$final
                state$currSeenLimit     <- 0
                state$currNotSeenLimit  <- 0
                state$numberOfReversals <- 0
                state$lastSeen          <- NA
            } else {
                state$finished <- TRUE
            }
        } else if (is.na(state$secondStairResult)) {
            state$secondStairResult <- state$lastSeen
            state$finished <- TRUE
        } else {
            warning("FT.step: stepping FT staircase when it has already terminated")
        }
    }

    return(list(state=state, resp=opiResp))
}#FT.step()

################################################################################
# Return TRUE if FT should stop, FALSE otherwise
#
# Input parameters
#   State list as returned by FT.start/step
# Returns 
#   TRUE or FALSE
################################################################################
FT.stop <- function(state) { return(state$finished) }

################################################################################
# Given a state, return an estimate of threshold and reason...
#
# Input parameters
#   State list as returned by FT.start/step
# Returns a list of
#   final      = final dB threshold
#   stopReason = reason staircase finished
################################################################################
FT.final.details <- function(state) {
    if (state$currSeenLimit == 2) {
        final <- state$maxStimulus
        stopReason <- "Max"
    } else if (state$currNotSeenLimit == 2) {
        stopReason <- "Min"
        final <- state$minStimulus
    } else {
        stopReason <- "Reversals"
        final <- state$lastSeen
    }
    
    return (list(
        final=final,                  # The threshold estimate in dB
        first=state$firstStairResult, # The threshold estimate in dB
        stopReason=stopReason,        # Reason for terminating staircase
        np=state$numPresentations
    ))  
}

################################################################################
# Given a state, return an estimate of threshold
#
# Input parameters
#   State list as returned by FT.start/step
# Returns 
#   final dB threshold
################################################################################
FT.final <- function(state) {
    return (FT.final.details(state)[["final"]])
}

########
# TEst

### require(OPI)
### chooseOpi("SimHenson")
### #chooseOpi("SimYes")
### #chooseOpi("SimNo")
### opiInitialize()
### makeStim <- function(db, n) { 
###     s <- list(x=9, y=9, level=dbTocd(db), size=0.43, color="white",
###              duration=1200, responseWindow=500)
###     class(s) <- "opiStaticStimulus"
### 
###     return(s)
### }
### 
### res <- lapply(0:40, function(tt) {
###     lapply(1:1000, function(i) {
###         s <- FT.start(makeStim=makeStim, tt=tt, fpr=0.15, fnr=0.3)
###         s <- FT.step(s)
###         while(!FT.stop(s$state)) {
###             s <- FT.step(s$state)
###         }
###         FT.final(s$state)
###     })
### })
