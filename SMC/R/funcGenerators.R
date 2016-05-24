
###  $Id: funcGenerators.R,v 1.2 2008/02/04 19:57:56 goswami Exp $
###  
###  File:    funcGenerators.R
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

## The so-called Markovian regime-switching model. Note the following
## function is common to the example sections of the help pages for
## particleFilter, auxiliaryParticleFilter and sequentialMonteCarlo,
## and so it returns more than what's required to run the example(s)
## in this help page.
MarkovSwitchingFuncGenerator <-
    function  (seed = -975313579)
{
    set.seed(seed)
    nTrain <- 100
    p00    <- 0.9
    p11    <- 0.95
    mu0    <- -1.5
    mu1    <- -mu0
    sigma  <- 0.75

    generateNextStreamsFunc <-
        function (currentPeriod,
                  lag1Streams,
                  lag1LogWeights,
                  streamIndices,
                  streamReps,
                  startingStreams = NULL,
                  ...)
        {
            if (currentPeriod == 1) {
                if (is.null(startingStreams)) {
                    return(as.matrix(rep(mu0, length(streamIndices))))
                }
                if ((nStartingStreams <- nrow(startingStreams)) !=
                    (nAuxiliaryVarIndices <- length(streamIndices))) {
                    msg <- paste('The startingStreams should be a matrix with',
                                 nAuxiliaryVarIndices, 'rows:')
                    fatal(msg, nStartingStreams)
                }
                return(startingStreams)               
            }

            ## currentPeriod > 1 case:
            lag1Streams <- lag1Streams[streamIndices, , drop = FALSE]
            ret         <- numeric(nrow(lag1Streams))

            filter0 <- (abs(lag1Streams[ , 1] - mu0) <= 0)
            if ((nn0 <- sum(filter0)) > 0) {
                ret[filter0] <- ifelse(runif(nn0) <= p00, mu0, mu1)
            }
            filter1 <- (abs(lag1Streams[ , 1] - mu1) <= 0)
            if ((nn1 <- sum(filter1)) > 0) {
                ret[filter1] <- ifelse(runif(nn1) <= p11, mu1, mu0)
            }
            as.matrix(ret)
        }
    
    ## generate the data
    mu     <- numeric(nTrain)
    yy     <- numeric(nTrain)
    tt     <- 1
    mu[tt] <- mu0
    yy[tt] <- rnorm(1, mu[tt], sigma)
    for (tt in seq.int(2, nTrain)) {
        if (abs(mu[tt - 1] - mu0) <= 0.0) {                
            mu[tt] <- ifelse((runif(1) <= p00), mu0, mu1)
        }
        else if (abs(mu[tt - 1] - mu1) <= 0.0) {
            mu[tt] <- ifelse((runif(1) <= p11), mu1, mu0)
        }
        else
            stop('should not reach here')
        
        yy[tt] <- rnorm(1, mu[tt], sigma)
    }
    
    logObsDensFunc <-
        function (currentPeriod, currentStreams, ...)
        {
            dnorm(yy[currentPeriod], mean = currentStreams[ , 1], sd = sigma,
                  log = TRUE)
        }
    
    generateStreamRepsFunc <-
        function (currentPeriod,
                  lag1Streams,
                  lag1LogWeights,
                  streamIndices,
                  ...)
        {
            if (currentPeriod == 1) return(rep(mu0, length(streamIndices)))
            c(sapply(lag1Streams,
                     function (ss)
                 {
                     if (abs(ss - mu0) <= 0.0)
                         return((1 - p00) * mu1 + p00 * mu0)
                     else if (abs(ss - mu1) <= 0.0)
                         return((1 - p11) * mu0 + p11 * mu1)
                     stop('should not reach here')                   
                 }))
        }

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
                    currentStreams <- as.matrix(rep(mu0, nStreamsToGenerate))
                } else {
                    if ((nStartingStreams <- nrow(startingStreams)) !=
                        nStreamsToGenerate) {
                        msg <- paste('The startingStreams should be a matrix with',
                                     nStreamsToGenerate, 'rows:')
                        fatal(msg, nStartingStreams)
                    }
                    currentStreams <- startingStreams[seq_len(nStreamsToGenerate), ]
                }
                currentLogDens <- dnorm(yy[currentPeriod],
                                        mean = currentStreams[ , 1],
                                        sd   = sigma,
                                        log  = TRUE)
                return(list(currentStreams    = currentStreams,
                            currentLogWeights = currentLogDens))
            }
            
            stopifnot(nStreamsToGenerate == nrow(lag1Streams))
            ret     <- numeric(nStreamsToGenerate)
            filter0 <- (abs(lag1Streams[ , 1] - mu0) <= 0)
            if ((nn0 <- sum(filter0)) > 0) {
                ret[filter0] <- ifelse(runif(nn0) <= p00, mu0, mu1)
            }
            filter1 <- (abs(lag1Streams[ , 1] - mu1) <= 0)
            if ((nn1 <- sum(filter1)) > 0) {
                ret[filter1] <- ifelse(runif(nn1) <= p11, mu1, mu0)
            }
            currentStreams <- as.matrix(ret)
            currentLogDens <- dnorm(yy[currentPeriod],
                                    mean = currentStreams[ , 1],
                                    sd   = sigma,
                                    log  = TRUE)
            return(list(currentStreams    = currentStreams,
                        currentLogWeights = currentLogDens))            
        }
    
    list(yy                      = as.matrix(yy),
         mu                      = as.matrix(mu),
         generateNextStreamsFunc = generateNextStreamsFunc,
         logObsDensFunc          = logObsDensFunc,
         generateStreamRepsFunc  = generateStreamRepsFunc,
         propagateFunc           = propagateFunc)    
}
