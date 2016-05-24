#
# Full Threshold (FT) algorithm for a single location.
# FT begins with a 4-2dB staircase at level est. If the final estimate 
# (last seen) is more than 4dB away from est, a second 4-2 staircase is
# completed beginning at the estimate returned from the first.
# There is special handling of the brightest and dimmest stimuli: see the 
# comments in doStair() below.
# Based on Andrew Turpin's implementation in his Barramundi Simulator.
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
#         Jonathan Denniss (jdenniss@unimelb.edu.au)
# Date: June 2012
#
# Copyright 2012 Andrew Turpin and Jonathan Denniss
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
#   verbose       True if you want each presentation printed
#   makeStim      A helper function to create the required
#                 OPI data type for passing to opiPresent
#   ...           Parameters for opiPresent
# Returns a list containing
#   npres    Total number of presentations
#   respSeq  Response sequence stored as a list of (seen,dB) pairs
#   first    First staircase estimate in dB
#   final    Final threshold estimate in dB
################################################################################
FT <- function(est=25, instRange=c(0,40), verbose=FALSE, makeStim, ...) {
    #
    # Do a single 4-2 staircase beginning at startEstimate
    # and stopping after 2 reversals, or if min/max not/seen twice.
    # Return the reason for stopping, the last seen (or min/max) 
    # and responseSequence
    #
    doStair <- function(startEstimate) {
        numRevs       <- 0    # number of reversals
        numPres       <- 0    # number of presentations
        numSeenMax    <- 0    # number of times seen instRange[2]
        numNotSeenMin <- 0    # number of times not seen instRange[1]
        lastSeen      <- NA   # last stimulus seen
        responseSeq   <- NULL # a list of (seen/not, db value) pairs

        currentEst <- startEstimate
        while (numRevs < 2 && numNotSeenMin < 2 && numSeenMax < 2) { 
            opiResp <- opiPresent(stim=makeStim(currentEst, numPres), nextStim=NULL, ...)
            while (!is.null(opiResp$err))
                opiResp <- opiPresent(stim=makeStim(currentEst, numPres), nextStim=NULL, ...)
            resp <- opiResp$seen
            numPres <- numPres + 1
            
            if (verbose) {
                cat(sprintf("Presentation %2d: ", numPres))
                cat(sprintf("dB= %2d repsonse=%s\n", currentEst, resp))
            }

            responseSeq <- c(responseSeq, list(c(seen=resp, db=currentEst)))

            if (resp)
                lastSeen <- currentEst

            if (currentEst == instRange[1] && !resp)
                numNotSeenMin <- numNotSeenMin + 1

            if (currentEst == instRange[2] && resp)
                numSeenMax <- numSeenMax + 1

            if (numPres > 1 && resp != responseSeq[[numPres-1]]["seen"])
                numRevs <-numRevs + 1 

            delta <- ifelse(numRevs == 0, 4, 2) * ifelse(resp, +1, -1)
            currentEst <- min(instRange[2], 
                            max(instRange[1], currentEst + delta))
        }

        if (numSeenMax == 2) {
            final <- instRange[2]
            stopReason <- "Max"
        } else if (numNotSeenMin == 2) {
            stopReason <- "Min"
            final <- instRange[1]
        } else {
            stopReason <- "Reversals"
            final <- lastSeen
        }

        return (list(
            stopReason=stopReason,     # Reason for terminating staircase
            final=final,               # The threshold estimate in dB
            responseSeq=responseSeq    # A list of (seen, db) pairs
        ))  
    }# doStair()

        #
        # do the first 4-2 stair and store the result
        #
    first <- doStair(est)
    final           <- first$final
    fullResponseSeq <- first$responseSeq

    if (first$stopReason == "Reversals" && abs(est - first$final) > 4) {
        second <- doStair(first$final)
        fullResponseSeq <- c(first$responseSeq, second$responseSeq)
        final <- second$final
    } 

    return(list(
        npres=length(fullResponseSeq),  # number of presentations
        respSeq=fullResponseSeq,        # reposnse sequence (list of pairs)
        first=first$final,              # estimate from first staircase
        final=final                     # final threshold estimate
    ))
}#FT()
