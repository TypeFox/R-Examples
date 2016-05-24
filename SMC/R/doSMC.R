
###  $Id: doSMC.R,v 1.21 2008/02/04 22:28:27 goswami Exp $
###  
###  File:    doSMC.R
###  Package: SMC
###  
###  Copyright (C) 2006-present Gopi Goswami
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 of the License, or
###  (at your option) any later version.
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  For a copy of the GNU General Public License please write to the
###  Free Software Foundation, Inc.
###  59 Temple Place, Suite 330.
###  Boston, MA  02111-1307 USA.
###
###  For bugs in the code please contact:
###  <goswami@stat.harvard.edu>
###
###  SYNOPSIS
###
###
###
###  DESCRIPTION
###
###
###

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### The following is the main workhorse function.

SMCMain <-
    function (nStreams,
              nPeriods,
              dimPerPeriod,
              propagateFunc, 
              resampCriterionFunc = NULL,
              resampFunc          = NULL,
              summaryFunc         = NULL,
              nMHSteps            = 0,
              MHUpdateFunc        = NULL,
              nStreamsPreResamp   = NULL, 
              returnStreams       = FALSE,
              returnLogWeights    = FALSE,
              verboseLevel        = 0,
              ...)   
{    
    ptm <- proc.time( )
    ## BEGIN: Error checks
    nStreams     <- .check.numericWithLLim(nStreams, 1)
    nPeriods     <- .check.numericWithLLim(nPeriods, 1)
    dimPerPeriod <- .check.numericWithLLim(dimPerPeriod, 1)

    propagateFunc <- .check.propagateFunc(propagateFunc)
    
    ## NOTE: The following functions could be passed as NULL, they
    ## have their builtin counterparts in C.
    resampCriterionFunc <- .check.resampCriterionFunc(resampCriterionFunc)
    resampFunc          <- .check.resampFunc(resampFunc)
    summaryFunc         <- .check.summaryFunc(summaryFunc)
    dimSummPerPeriod    <- .get.dimSummPerPeriod(nStreams     = nStreams,
                                                 dimPerPeriod = dimPerPeriod,
                                                 summaryFunc  = summaryFunc,
                                                 ...)
    
    nMHSteps <- .check.numericWithLLim(nMHSteps, 0)
    if (nMHSteps > 0)
        MHUpdateFunc <- .check.MHUpdateFunc(MHUpdateFunc)
    else
        MHUpdateFunc <- NULL
    
    nStreamsPreResamp <- .check.nStreamsPreResamp(nStreamsPreResamp, nStreams)
    returnStreams     <- .check.logical(returnStreams)
    returnLogWeights  <- .check.logical(returnLogWeights)
    verboseLevel      <- .check.numericWithLLim(verboseLevel, NA)  
    procTimeFunc      <- as.function(proc.time)             
    procTimeFuncEnv   <- new.env( )
    doCallFunc        <- as.function(doCall)             
    doCallFuncEnv     <- new.env( )
    dotsList          <- list(...)
    argsList          <- collectVarnames(ls( ))
    ## E N D: Error checks

    if (argsList$verboseLevel >= 3) {
        cat('The processed arguments:\n')
        print(argsList)
    }
    if (argsList$verboseLevel >= 1) cat('\nBEGIN: SMC\n')        
    res <- .Call('SMCMainC', argsList)
    if (argsList$verboseLevel >= 1) cat('E N D: SMC\n')

    res$nStreams          <- nStreams
    res$nPeriods          <- nPeriods
    res$dimPerPeriod      <- dimPerPeriod
    res$nStreamsPreResamp <- nStreamsPreResamp
    res$nMHSteps          <- nMHSteps
    res                   <- c(res, list(summaryFunc = summaryFunc))

    if (!any(is.na(ptm))) res$time <- proc.time( ) - ptm
    class(res) <- 'SMC'
    res
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequentialMonteCarlo <-
    function (nStreams,
              nPeriods,
              dimPerPeriod,
              propagateFunc, 
              resampCriterionFunc = NULL,
              resampFunc          = NULL,
              summaryFunc         = NULL,
              nMHSteps            = 0,
              MHUpdateFunc        = NULL,
              nStreamsPreResamp   = NULL, 
              returnStreams       = FALSE,
              returnLogWeights    = FALSE,
              verboseLevel        = 0,
              ...)                          
{

    ret <- SMCMain(nStreams            = nStreams,
                   nPeriods            = nPeriods,
                   dimPerPeriod        = dimPerPeriod,
                   propagateFunc       = propagateFunc,
                   resampCriterionFunc = resampCriterionFunc,
                   resampFunc          = resampFunc,
                   summaryFunc         = summaryFunc,
                   nMHSteps            = nMHSteps,
                   MHUpdateFunc        = MHUpdateFunc,
                   nStreamsPreResamp   = nStreamsPreResamp,
                   returnStreams       = returnStreams,
                   returnLogWeights    = returnLogWeights,
                   verboseLevel        = verboseLevel,
                   ...)
    ret$propagateFunc <- propagateFunc
    ret$filterType    <- 'sequentialMonteCarlo'
    ret
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Some helper functions to be used by the rest of the exported
### functions in the package.

.get.streamsMat <-
    function (streams)
{
    if (is.vector(streams))
        return(as.matrix(streams))

    if (is.matrix(streams))
        return(t(streams))

    msg <- 'The provided streams should be either a vector or a matrix'
    fatal(msg, streams)
}

.get.probs <-
    function (logWeights)
{
    weights <- exp(logWeights - max(logWeights))
    weights / sum(weights)
}

.resampAlways <- 
    function (currentPeriod,
              currentStreams,
              currentLogWeights,
              ...)
{ TRUE }

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Some helper functions for auxiliaryParticleFilter

.propagateFuncGen.APF <-
    function (gSRF, gNSF, lODF)
{
    propagateFunc <-
        function (currentPeriod,
                  nStreamsToGenerate,
                  lag1Streams,
                  lag1LogWeights,
                  startingStreams = NULL,
                  ...)
        {
            if (currentPeriod == 1) {
                if (is.null(startingStreams)) {
                    currentStreams <-
                        gNSF(currentPeriod   = currentPeriod,
                             lag1Streams     = NULL,
                             lag1LogWeights  = NULL,
                             streamIndices   = seq_len(nStreamsToGenerate),
                             streamReps      = NULL,
                             startingStreams = startingStreams,
                             ...)
                }
                else {
                    if (!is.matrix(startingStreams) ||
                        ((nStartingStreams <- nrow(startingStreams)) <
                         nStreamsToGenerate)) {
                        msg <- paste('The startingStreams should be a matrix',
                                     'with >=', nStreamsToGenerate, 'rows:')
                        fatal(msg, nStartingStreams)
                    }
                    
                    currentStreams <- startingStreams[seq_len(nStreamsToGenerate), ]
                }
                
                currentLogDens <- lODF(currentPeriod, currentStreams, ...)
                return(list(currentStreams    = currentStreams,
                            currentLogWeights = currentLogDens))
            }

            ## currentPeriod > 1 case:
            ## generate the streamReps
            lag1Streams <- as.matrix(lag1Streams)
            nStreams    <- nrow(lag1Streams)
            mus         <- as.matrix(gSRF(currentPeriod  = currentPeriod,
                                          lag1Streams    = lag1Streams,
                                          lag1LogWeights = lag1LogWeights,
                                          streamIndices  = seq_len(nStreams),
                                          ...))

            ## sample the streamReps
            musLogDens    <- lODF(currentPeriod, mus, ...)
            probsKK       <- .get.probs(lag1LogWeights + musLogDens)
            sampKK        <- sample(seq_len(nStreams),
                                    size    = nStreamsToGenerate,
                                    replace = TRUE,
                                    prob    = probsKK)

            ## propagate the sampled streams, compute the logDens on
            ## these streams and return
            sampStreams     <- gNSF(currentPeriod   = currentPeriod,
                                    lag1Streams     = lag1Streams,
                                    lag1LogWeights  = lag1LogWeights,
                                    streamIndices   = sampKK,
                                    streamReps      = mus,
                                    startingStreams = startingStreams,
                                    ...)
            sampLogDens <- lODF(currentPeriod, sampStreams, ...)

            list(currentStreams    = sampStreams,
                 currentLogWeights = sampLogDens - musLogDens[sampKK])
        }

    propagateFunc
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auxiliaryParticleFilter <-
  function (nStreams,
            nPeriods,
            dimPerPeriod,
            generateStreamRepsFunc,
            generateNextStreamsFunc,
            logObsDensFunc,
            resampCriterionFunc = NULL,
            resampFunc          = NULL,
            summaryFunc         = NULL,
            nMHSteps            = 0,
            MHUpdateFunc        = NULL,
            nStreamsPreResamp   = NULL,
            returnStreams       = FALSE,
            returnLogWeights    = FALSE,
            verboseLevel        = 0,
            ...)                                 
{
    generateStreamRepsFunc  <-
        .check.generateStreamRepsFunc(generateStreamRepsFunc)
    generateNextStreamsFunc <-
        .check.generateNextStreamsFunc.APF(generateNextStreamsFunc)
    logObsDensFunc          <-
        .check.logObsDensFunc(logObsDensFunc)
    propagateFunc           <-
        .propagateFuncGen.APF(gSRF = generateStreamRepsFunc,
                              gNSF = generateNextStreamsFunc,
                              lODF = logObsDensFunc)

    ## The following is the default setting specific to APF
    if (is.null(nStreamsPreResamp))
        nStreamsPreResamp <- as.integer(nStreams * 1.25)
    
    ## In APF, we need to always resample since:
    ## nStreamsPreResamp > nStreams
    resampCriterionFunc <- .resampAlways
    
    ret <- SMCMain(nStreams            = nStreams,
                   nPeriods            = nPeriods,
                   dimPerPeriod        = dimPerPeriod,
                   propagateFunc       = propagateFunc,
                   resampCriterionFunc = resampCriterionFunc, 
                   resampFunc          = resampFunc,
                   summaryFunc         = summaryFunc,
                   nMHSteps            = nMHSteps,
                   MHUpdateFunc        = MHUpdateFunc,
                   nStreamsPreResamp   = nStreamsPreResamp,
                   returnStreams       = returnStreams,
                   returnLogWeights    = returnLogWeights,
                   verboseLevel        = verboseLevel,
                   ...)
    ret$generateStreamRepsFunc  <- generateNextStreamsFunc
    ret$generateNextStreamsFunc <- generateNextStreamsFunc
    ret$logObsDensFunc          <- logObsDensFunc
    ret$filterType              <- 'auxiliaryParticleFilter'
    ret    
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Some helper functions for particleFilter

.propagateFuncGen.PF <-
    function (gNSF, lODF)
{
    propagateFunc <-
        function (currentPeriod,
                  nStreamsToGenerate,
                  lag1Streams,
                  lag1LogWeights,
                  startingStreams = NULL,
                  ...)
        {
            if (currentPeriod == 1) {
                if (is.null(startingStreams)) {
                    currentStreams <-
                        gNSF(currentPeriod   = currentPeriod,
                             lag1Streams     = NULL,
                             lag1LogWeights  = NULL,
                             streamIndices   = seq_len(nStreamsToGenerate),
                             startingStreams = startingStreams,
                             ...)
                }
                else {
                    if (!is.matrix(startingStreams) ||
                        ((nStartingStreams <- nrow(startingStreams)) <
                         nStreamsToGenerate)) {
                        msg <- paste('The startingStreams should be a matrix',
                                     'with >=', nStreamsToGenerate, 'rows:')
                        fatal(msg, nStartingStreams)
                    }
                    
                    currentStreams <- startingStreams[seq_len(nStreamsToGenerate), ]
                }
                
                currentLogDens <- lODF(currentPeriod, currentStreams, ...)
                return(list(currentStreams    = currentStreams,
                            currentLogWeights = currentLogDens))
            }

            if ((nLag1Streams <- nrow(lag1Streams)) != nStreamsToGenerate) {
                msg <- paste('Given nStreamsToGenerate (= ', nStreamsToGenerate,
                             ') != nLag1Streams (= ', nLag1Streams, ') for ',
                             'currentPeriod = ', currentPeriod, '. I don\'t know ',
                             'how to handle this case in the context of particle ',
                             'filtering. Please implement your idea through the ',
                             'function "generateNextStreamsFunc" and rerun.')
                fatal(msg)
            }
            currentStreams <- gNSF(currentPeriod   = currentPeriod,
                                   lag1Streams     = lag1Streams,
                                   lag1LogWeights  = lag1LogWeights,
                                   streamIndices   = seq_len(nStreamsToGenerate),
                                   startingStreams = startingStreams,
                                   ...)
            currentLogDens <- lODF(currentPeriod, currentStreams, ...)

            list(currentStreams    = currentStreams,
                 currentLogWeights = lag1LogWeights + currentLogDens)
        }

    propagateFunc
}

.MHUpdateFuncGen.PF <-
    function (gNSF, lODF)
{
    MHUpdateFunc <-
        function (currentPeriod,
                  nMHSteps,
                  currentStreams,
                  lag1Streams,
                  lag1LogWeights,
                  ...)
        {
            stopifnot(currentPeriod >= 1)
            currentStreams <- as.matrix(currentStreams)
            nStreams       <- nrow(currentStreams)
            rr             <- seq_len(nStreams)
            prob           <- .get.probs(lag1LogWeights)
            currentLogDens <- lODF(currentPeriod, currentStreams, ...)
            accepted       <- numeric(nStreams)
            
            for (ii in seq_len(nMHSteps)) {
                ## Sample the lag1Streams to be used to generate
                ## proposals
                streamIndices    <- sample(rr,
                                           size    = nStreams,
                                           replace = TRUE,
                                           prob    = prob)

                ## Generate the proposals using the chosen lag1Streams
                proposedStreams <- gNSF(currentPeriod   = currentPeriod,
                                        lag1Streams     = lag1Streams,
                                        lag1LogWeights  = lag1LogWeights,
                                        streamIndices   = streamIndices,
                                        startingStreams = NULL,
                                        ...)
                proposedLogDens  <- lODF(currentPeriod, proposedStreams, ...)

                ## The acceptance / rejection step
                decision                   <- (log(runif(nStreams)) <=
                                               (proposedLogDens - currentLogDens))
                accepted                   <- accepted + decision
                currentStreams[decision, ] <- proposedStreams[decision, ]
                currentLogDens[decision]   <- proposedLogDens[decision]
            }

            list(currentStreams  = currentStreams,
                 acceptanceRates = accepted / nMHSteps)
        }

    MHUpdateFunc
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

particleFilter <-
  function (nStreams,
            nPeriods,
            dimPerPeriod,
            generateNextStreamsFunc,
            logObsDensFunc,
            resampCriterionFunc = NULL,
            resampFunc          = NULL,
            summaryFunc         = NULL,
            nMHSteps            = 0,
            MHUpdateFunc        = NULL,
            nStreamsPreResamp   = NULL,
            returnStreams       = FALSE,
            returnLogWeights    = FALSE,
            verboseLevel        = 0,
            ...)                                
{
    generateNextStreamsFunc <-
        .check.generateNextStreamsFunc.PF(generateNextStreamsFunc)
    logObsDensFunc          <-
        .check.logObsDensFunc(logObsDensFunc)
    propagateFunc           <-
        .propagateFuncGen.PF(gNSF = generateNextStreamsFunc,
                             lODF = logObsDensFunc)

    if (is.null(MHUpdateFunc))
        MHUpdateFunc <- .MHUpdateFuncGen.PF(gNSF = generateNextStreamsFunc,
                                            lODF = logObsDensFunc)

    ## If nStreamsPreResamp > nStreams then always resamp
    nStreams            <- .check.numericWithLLim(nStreams, 1)
    nStreamsPreResamp   <- .check.nStreamsPreResamp(nStreamsPreResamp, nStreams)
    if (nStreamsPreResamp > nStreams)
        resampCriterionFunc <- .resampAlways

    ret <- SMCMain(nStreams            = nStreams,
                   nPeriods            = nPeriods,
                   dimPerPeriod        = dimPerPeriod,
                   propagateFunc       = propagateFunc,
                   resampCriterionFunc = resampCriterionFunc, 
                   resampFunc          = resampFunc,
                   summaryFunc         = summaryFunc,
                   nMHSteps            = nMHSteps,
                   MHUpdateFunc        = MHUpdateFunc,
                   nStreamsPreResamp   = nStreamsPreResamp,
                   returnStreams       = returnStreams,
                   returnLogWeights    = returnLogWeights,
                   verboseLevel        = verboseLevel,
                   ...)
    ret$generateNextStreamsFunc <- generateNextStreamsFunc
    ret$logObsDensFunc          <- logObsDensFunc
    ret$filterType              <- 'particleFilter'
    ret        
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.SMC <-
    function (x, ...)
{
    varnames <- c('nStreams', 'nStreamsPreResamp', 'nPeriods', 'dimPerPeriod')
    if (x$nMHSteps > 0)
        varnames <- c(varnames, 'nMHSteps')
    tmp           <- t(data.frame(x[varnames]))
    rownames(tmp) <- paste(rownames(tmp), ':', sep = '')
    colnames(tmp) <- ''
    cat('\nThe settings for this', x$filterType, 'run:\n')
    print(tmp)
    tmp <- x$draws$propUniqueStreamIds
    cat('\nThe summary of the resampling proportions:\n',
        '[Note: resampling was done for ', sum(!is.na(tmp)), ' periods]\n\n',
        sep = '')
    print(summary(tmp))
    cat('\n')
}

