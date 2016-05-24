#
# 4-2 staircase with 2 reversals.
# See comments in full_threshold.r.
# This does not initiate a second staircase.
# This returns averge of last two presentations as threshold.
#
# Includes
#     fourTwo.start    # initialise list state
#     fourTwo.step     # take state, present stim, update and return state
#     fourTwo.stop     # boolean - true if state is finished
#     fourTwo.final    # return final estimate from state
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: December 2014
#
# Copyright 2014 Andrew Turpin
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
#   verbose       If true, print data for each presentation
#   makeStim      A helper function to create the required
#                 OPI data type for passing to opiPresent
#   ...           Parameters for opiPresent
# Returns a list containing
#   npres    Total number of presentations
#   respSeq  Response sequence stored as a list of (seen,dB) pairs
#   final    Final threshold estimate in dB
################################################################################
fourTwo.start <- function(est=25, instRange=c(0,40), verbose=FALSE, makeStim, ...) {
    if (est < instRange[1] || est > instRange[2])
        stop("fourTwo.start: est must be in the range of instRange")

    return(list(name="fourTwo",
        startingEstimate=est,
        currentLevel=est,
        minStimulus=instRange[1],
        maxStimulus=instRange[2],
        makeStim=makeStim,
        lastSeen=NA,
        lastResponse=NA,
        stairResult=NA,
        finished="Not",                     # "Not", or one of "Max", "Min", "Rev"
        verbose=verbose,
        numberOfReversals=0,
        currSeenLimit=0,                    # number of times maxStimulus seen
        currNotSeenLimit=0,                 # number of times minStimulus not seen
        numPresentations=0,                 # number of presentations so far
        stimuli=NULL,                       # vector of stims shown
        responses=NULL,                     # vector of responses (1 seen, 0 not)
        responseTimes=NULL,                 # vector of response times
        opiParams=list(...)                 # the extra params
    ))
}# fourTwo.start()

################################################################################
# Present state$currentLevel and update state after response.
#
# Input parameters
#   State list as returned by fourTwo.start
#   nextStim - suitable for passing to opiPresent, can be NULL
# Returns a list containing
#   state = updated state list
#   resp  = response to the presentation (as returned by opiPresent)
#
# Follows fourTwo algorithm as documented in full_threshold.r.
# Note 
#   1) opiPresent called infinitely until no error
################################################################################
fourTwo.step <- function(state, nextStim=NULL) {
    if (state$finished != "Not")
        warning("fourTwo.step: stepping fourTwo staircase when it has already terminated")

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

    if (state$verbose) {
        cat(sprintf("Presentation %2d: ", state$numPresentations))
        cat(sprintf("dB= %2d repsonse=%s\n", state$currentLevel, opiResp$seen))
    }

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
    if (state$numberOfReversals >= 2) {
        state$finished <- "Rev"
        state$stairResult <- mean(tail(state$stimuli, 2)) # mean of last two
    } else if (state$currNotSeenLimit >= 2) {
        state$finished <- "Min"
        state$stairResult <- state$minStimulus
    } else if (state$currSeenLimit >= 2) {
        state$finished <- "Max"
        state$stairResult <- state$maxStimulus
    } else {
            # update stimulus for next presentation
        delta <- ifelse(state$numberOfReversals == 0, 4, 2) * ifelse(opiResp$seen, +1, -1)
        state$currentLevel <- min(state$maxStimulus, max(state$minStimulus, state$currentLevel + delta))
    } 

    return(list(state=state, resp=opiResp))
}#fourTwo.step()

################################################################################
# Return TRUE if fourTwo should stop, FALSE otherwise
#
# Input parameters
#   State list as returned by fourTwo.start/step
# Returns 
#   TRUE or FALSE
################################################################################
fourTwo.stop <- function(state) { return(state$finished != "Not") }

################################################################################
# Given a state, return an estimate of threshold
#
# Input parameters
#   State list as returned by fourTwo.start/step
# Returns 
#   final dB threshold
################################################################################
fourTwo.final <- function(state) {
    if (state$finished != "Not")
        return (state$stairResult)
    else {
        warning("fourTwo.step: asking for final result of unfinished staircase")
        return(NA)
    }
}

########
# Test
###
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
###         s <- fourTwo.start(makeStim=makeStim, tt=tt, fpr=0.15, fnr=0.3)
###         s <- fourTwo.step(s)
###         while(!fourTwo.stop(s$state)) {
###             s <- fourTwo.step(s$state)
###         }
###         fourTwo.final(s$state)
###     })
### })
