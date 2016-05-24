#
# ZEST algorithm that maintains state in a list - good for interleaved.
# Includes
#     ZEST          # just for a single location
#     ZEST.start    # initialise list state
#     ZEST.step     # take state, present stim, update and return state
#     ZEST.stop     # boolean - true if state is finished
#     ZEST.final    # return final estimate from state
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: August 2012
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

##########################
# little helper functions
##########################
ZEST.stdev <- function(state) { sqrt(sum(state$pdf*state$domain*state$domain) - sum(state$pdf * state$domain)^2) }

ZEST.entropy <- function(state) {
    z <- which(state$pdf > 0)
    return(-sum(state$pdf[z] * log2(state$pdf[z])))
}


##############################################
# Initialise state list
#   domain        List of dB values over which pdf is kept
#   prior         Probability distribution over domain.
#   likelihood    matrix where likelihood[s,t] is likelihood of seeing s given t is true thresh (Pr(s|t)
#                 where s and t are indexs into domain
#   stopType      N | S | H
#   stopValue     Value for num prs (N), stdev (S) of Entropy (H)
#   minStimulus   Lowest value to present
#   maxStimulus   Highest value to present
#   minNotSeenLimit Will terminate if minLimit value not seen this many times
#   maxSeenLimit    Will terminate if maxLimit value seen this many times
#   maxPresentations Maximum number of presentations
#   makeStim      A helper function to create the required
#                 OPI data type for passing to opiPresent
#   stimChoice    "mean", "median", "mode"
#   ...           Parameters for opiPresent
##############################################
ZEST.start <- function(domain=0:40, prior=rep(1/length(domain),length(domain)), 
            likelihood=sapply(domain, function(tt) { 0.03 + (1-0.03-0.03)*(1-pnorm(domain, tt, 1)) }),
            stopType="S",
            stopValue=1.5,
            minStimulus=head(domain, 1),
            maxStimulus=tail(domain, 1),
            maxSeenLimit=2,
            minNotSeenLimit=2,
            maxPresentations=100,
            makeStim, 
            stimChoice="mean",
            ...) {
    ##########################
    # Validate params
    ##########################
    if (!is.element(stopType, c("S", "H", "N")))
        stop("ZEST.start: stopType must be one of 'S', 'N', or 'H'")
    if (nrow(likelihood) != length(domain))
        stop(paste("ZEST.start: not enough rows in likelihood. Expect",length(domain)))
    if (ncol(likelihood) != length(domain))
        stop(paste("ZEST.start: not enough cols in likelihood. Expect",length(domain)))

    if (!is.element(minStimulus, domain))
        warning(paste("ZEST.start: you specified minStimulus=",minStimulus,"but it is not in domain."))
    if (!is.element(maxStimulus, domain))
        warning(paste("ZEST.start: you specified maxStimulus=",maxStimulus,"but it is not in domain."))

    pdf <- prior/sum(prior)

    return(list(name="ZEST",
                domain=domain, 
                pdf=pdf,
                likelihood=likelihood, 
                stopType=stopType,
                stopValue=stopValue,
                minStimulus=minStimulus,
                maxStimulus=maxStimulus,
                maxSeenLimit=maxSeenLimit,
                minNotSeenLimit=minNotSeenLimit,
                maxPresentations=maxPresentations,
                makeStim=makeStim,
                stimChoice=stimChoice,
                currSeenLimit=0,                    # number of times maxStimulus seen
                currNotSeenLimit=0,                 # number of times minStimulus not seen
                numPresentations=0,                 # number of presentations so far
                stimuli=NULL,                       # vector of stims shown
                responses=NULL,                     # vector of responses (1 seen, 0 not)
                responseTimes=NULL,                 # vector of response times
                opiParams=list(...)                 # the extra params
            ))
}# ZEST.start

################################################################################
# Present current stim and update state after response.
#
# Input parameters
#   State list as returned by ZEST.start
#   nextStim - suitable for passing to opiPresent, can be NULL
# Returns a list containing
#   state = updated state list
#   resp  = response to the presentation (as returned by opiPresent)
#
# Note 
#   1) stims are rounded to nearest domain entry 
#   2) opiPresent called infinitely until no error
################################################################################
ZEST.step <- function(state, nextStim=NULL) {

    if (state$stimChoice == "mean") {
        stimIndex <- which.min(abs(state$domain - sum(state$pdf * state$domain)))
    } else if (state$stimChoice == "mode") {
        stimIndex <- which.max(state$pdf)
    } else if (state$stimChoice == "median") {
        stimIndex <- which.min(abs(cumsum(state$pdf) - 0.5))
    } else {
        stop(paste("ZEST.step: stimChoice = ",state$stimChoice," not implemented."))
    }
    stim <- state$domain[stimIndex]
    stim <- max(stim, state$minStimulus)   # check not outside [minStimulus,maxStimulus]
    stim <- min(stim, state$maxStimulus)

    if (is.null(state$opiParams))
        params <- list(stim=state$makeStim(stim, state$numPresentations), nextStim=nextStim)
    else
        params <- c(list(stim=state$makeStim(stim, state$numPresentations), nextStim=nextStim), state$opiParams)
    opiResp <- do.call(opiPresent, params)
    while (!is.null(opiResp$err))
        opiResp <- do.call(opiPresent, params)
    state$stimuli          <- c(state$stimuli, stim)
    state$responses        <- c(state$responses, opiResp$seen)
    state$responseTimes    <- c(state$responseTimes, opiResp$time)
    state$numPresentations <- state$numPresentations + 1
    
    if(opiResp$seen) { 
        if (stim == state$maxStimulus) state$currSeenLimit <- state$currSeenLimit + 1
        state$pdf <- state$pdf * state$likelihood[stimIndex, ]
    } else {
        if (stim == state$minStimulus) state$currNotSeenLimit <- state$currNotSeenLimit + 1
        state$pdf <- state$pdf * (1 - state$likelihood[stimIndex, ])
    }
    state$pdf <- state$pdf/sum(state$pdf)

    return(list(state=state, resp=opiResp))
}#ZEST.step()

################################################################################
# Return TRUE if ZEST should stop, FALSE otherwise
#
# Input parameters
#   State list as returned by ZEST.start/step
# Returns 
#   TRUE or FALSE
################################################################################
ZEST.stop <- function(state) {
    keepGoing <- (
        (state$numPresentations < state$maxPresentations) &&
        (state$currNotSeenLimit < state$minNotSeenLimit) &&
        (state$currSeenLimit    < state$maxSeenLimit) &&
        (
           ((state$stopType == "S") && (ZEST.stdev(state) > state$stopValue))
        || ((state$stopType == "H") && (ZEST.entropy(state) > state$stopValue))
        || ((state$stopType == "N") && (state$numPresentations < state$stopValue))
        )
    )
    return (!keepGoing)
}#ZEST.stop

################################################################################
# Given a state, return an estimate of threshold
#
# Input parameters
#   State list as returned by ZEST.start/step
# Returns 
#   Mean   of pdf if state$stimChoice == "mean"
#   Mode   of pdf if state$stimChoice == "mode"
#   Median of pdf if state$stimChoice == "median"
################################################################################
ZEST.final <- function(state) {
    if (state$stimChoice == "mean") {
        final <- sum(state$pdf*state$domain)
    } else if (state$stimChoice == "mode") {
        final <- state$domain[which.max(state$pdf)]
    } else if (state$stimChoice == "median") {
        final <- state$domain[which.min(abs(cumsum(state$pdf) - 0.5))]
    } 

    return(final)
}#ZEST.final

################################################################################
# ZEST for a single location. 
#
# Input parameters
#   domain        List of dB values over which pdf is kept
#   prior         Probability distribution over domain.
#   likelihood    matrix where likelihood[s,t] is likelihood of seeing s given t is true thresh (Pr(s|t)
#                 where s and t are indexs into domain
#   stopType      N | S | H
#   stopValue     Value for num prs (N), stdev (S) of Entropy (H)
#   minNotSeenLimit Will terminate if lowest domain value not seen this many times
#   maxSeenLimit    Will terminate if highest domain value seen this many times
#   maxPresentations Maximum number of presentations
#   verbose       1 if you want pdfs returned, 2 is 1+print, 0 for none
#   makeStim      A helper function to create the required
#                 OPI data type for passing to opiPresent
#   stimChoice    "mean", "median", "mode"
#   ...           Parameters for opiPresent
# Returns a list containing
#   npres    Total number of presentations
#   respSeq  Response sequence stored as a matrix: row 1 = dB, row 2 = response 1/0
#   pdfs     Sequence of pdfs used (if verbose)
#
# Note 
#   1) stims are rounded to nearest domain entry 
#   2) opiPresent called infinitely until no error
################################################################################
ZEST <- function(domain=0:40, prior=rep(1/length(domain),length(domain)), 
            likelihood=sapply(domain, function(tt) { 0.03 + (1-0.03-0.03)*(1-pnorm(domain, tt, 1)) }),
            stopType="S",
            stopValue=1.5,
            minStimulus=head(domain,1),
            maxStimulus=tail(domain,1),
            maxSeenLimit=2,
            minNotSeenLimit=2,
            maxPresentations=100,
            verbose=0, makeStim, 
            stimChoice="mean",
            ...) {
    state <- ZEST.start(domain, prior, likelihood, stopType, stopValue, 
                        minStimulus, maxStimulus, 
                        maxSeenLimit,minNotSeenLimit,maxPresentations,
                        makeStim,stimChoice, ...)

    pdfs <- NULL
    while(!ZEST.stop(state)) {
        r <- ZEST.step(state)
        state <- r$state
        if (verbose == 2) {
            cat(sprintf("Presentation %2d: ", state$numPresentations))
            cat(sprintf("stim= %5s repsonse=%s ", tail(state$stimuli,1), tail(state$responses,1)))
            cat(sprintf("stdev= %8.4g H= %8.4g\n", ZEST.stdev(state), ZEST.entropy(state)))
        }
        if (verbose > 0)
            pdfs <- c(pdfs, list(state$pdf))
    }

    return(list(
        npres=tail(state$numPresentations,1),        # number of presentations
        respSeq=mapply(c, state$stimuli, state$responses), # reposnse sequence (list of pairs)
        pdfs=pdfs,                                   # list of pdfs used (if verbose > 0)
        final=ZEST.final(state)                      # final threshold estimate
    ))
}#ZEST

############################################################
# Tests
############################################################
#require(OPI)
#chooseOpi("SimHenson")
##chooseOpi("SimYes")
#opiInitialize("C",6)
#
#makeStim <- function(db, n) { 
#         s <- list(x=9, y=9, level=dbTocd(db,10000/pi), size=0.43, 
#                  color="white",
#                  duration=200, responseWindow=1500)
#         class(s) <- "opiStaticStimulus"
#     
#         return(s)
#     }
#makeNextStim <- function(x,y) { 
#         s <- list(x=9, y=9, level=dbTocd(db,10000), size=0.43, color="white",
#                  duration=200, responseWindow=1500)
#         class(s) <- "opiStaticStimulus"
#     
#         return(s)
#     }
#
#state <- ZEST.start(domain=-5:45, maxStimulus=40, minStimulus=0, makeStim=makeStim, stopType="S", stopValue= 1.5, tt=0, fpr=0.30)
#while(!ZEST.stop(state)) {
#    r <- ZEST.step(state)
#    cat(sprintf("%2d %s\n",tail(r$state$stimuli,1), r$resp$seen))
#    state <- r$state
#}
#print(ZEST.final(state))
#
#print("########################################################################")
#
#makeStimHelper <- function(db,n, x, y) {
#    ff <- function(db, n) db+n
#
#    body(ff) <- substitute(
#        {s <- list(x=x, y=y, level=dbTocd(db,10000), size=0.43, color="white",
#                  duration=200, responseWindow=1500)
#         class(s) <- "opiStaticStimulus"
#         return(s)
#        }
#        , list(x=x,y=y)) 
#    return(ff)
#} 
#
#    # list of (x, y, true threshold) triples
#locations <- list(c(9,9,30), c(-9,-9,32), c(9,-9,31), c(-9,9,33))
#
#    # setup starting states for each location
#states <- lapply(locations, function(loc) {
#    ZEST.start(domain=-5:45,
#        makeStim=makeStimHelper(db,n,loc[1],loc[2]),
#        maxStimulus=40, minStimulus=0,         
#        stopType="S", stopValue= 1.5, tt=loc[3], fpr=0.03, fn=0.01)
#})
#
#    # loop through until all states are "stop"
#while(!all(st <- unlist(lapply(states, ZEST.stop)))) {
#    i <- sample(which(!st), 1)  # choose a random, unstopped state
#    r <- ZEST.step(states[[i]]) # step it
#    states[[i]] <- r$state      # update the states
#}
#
#finals <- lapply(states, ZEST.final)    # get final estimates of threshold
#for(i in 1:length(locations))
#    cat(sprintf("Location (%+2d,%+2d) has threshold %4.2f\n",locations[[i]][1], locations[[i]][2], finals[[i]]))
#    
#print("########################################################################")
#
#a <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stopType="H", stopValue=  3, verbose=0, tt=20, fpr=0.03))
#b <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stopType="S", stopValue=1.5, verbose=0, tt=20, fpr=0.03))
#c <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stopType="S", stopValue=2.0, verbose=0, tt=20, fpr=0.03))
#d <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stopType="N", stopValue= 50, verbose=0, tt=20, fpr=0.03))
#
#a <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stimChoice="mean", tt=20, fpr=0.03))
#b <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stimChoice="mode", tt=20, fpr=0.03))
#c <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stimChoice="median", tt=20, fpr=0.03))
#d <- sapply(1:100, function(i) ZEST(makeStim=makeStim, stimChoice="mean", tt=20, fpr=0.03))
#
#layout(matrix(1:2,1,2))
#boxplot(lapply(list(a,b,c,d), function(x) unlist(x["final",])))
#boxplot(lapply(list(a,b,c,d), function(x) unlist(x["npres",])))
#
#a <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=0,s=5), tt=30, fpr=0.03))
#b <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=10,s=5), tt=30, fpr=0.03))
#c <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=20,s=5), tt=30, fpr=0.03))
#d <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=30,s=5), tt=30, fpr=0.03))
#layout(matrix(1:2,1,2))
#boxplot(lapply(list(a,b,c,d), function(x) unlist(x["final",])))
#boxplot(lapply(list(a,b,c,d), function(x) unlist(x["npres",])))
#
#repp <- function(...) sapply(1:100, function(i) ZEST(makeStim=makeStim, ...))
#a <- repp(stopType="H", stopValue=  3, verbose=0, tt=30, fpr=0.03)
#b <- repp(stopType="S", stopValue=1.5, verbose=0, tt=30, fpr=0.03)
#c <- repp(stopType="S", stopValue=2.0, verbose=0, tt=30, fpr=0.03)
#d <- repp(stopType="N", stopValue= 50, verbose=0, tt=30, fpr=0.03)
#e <- repp(prior=dnorm(0:40,m=0,s=5), tt=30, fpr=0.03)
#f <- repp(prior=dnorm(0:40,m=10,s=5), tt=30, fpr=0.03)
#g <- repp(prior=dnorm(0:40,m=20,s=5), tt=30, fpr=0.03)
#h <- repp(prior=dnorm(0:40,m=30,s=5), tt=30, fpr=0.03)
#
#layout(matrix(1:2,1,2))
#boxplot(lapply(list(a,b,c,d,e,f,g,h), function(x) unlist(x["final",])))
#boxplot(lapply(list(a,b,c,d,e,f,g,h), function(x) unlist(x["npres",])))
#
#a <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=30,s=5), tt=00, fpr=0.03))
#b <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=30,s=5), tt=10, fpr=0.03))
#c <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=30,s=5), tt=20, fpr=0.03))
#d <- sapply(1:100, function(i) ZEST(makeStim=makeStim, prior=dnorm(0:40,m=30,s=5), tt=30, fpr=0.03))
#layout(matrix(1:2,1,2))
#boxplot(lapply(list(a,b,c,d), function(x) unlist(x["final",])))
#boxplot(lapply(list(a,b,c,d), function(x) unlist(x["npres",])))
#
#
