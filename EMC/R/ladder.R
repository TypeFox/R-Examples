

###  $Id: ladder.R,v 1.25 2008/02/03 04:18:52 goswami Exp $
###  
###  File:    ladder.R
###  Package: EMC
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


### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Function for the Parzen kernel.
### See Andrews, D, Econometrica, 1991, P821

kernelParzen <- function (xx)
{
    yy <- abs(xx)
    if (yy > 1.0) {
        return(0.0)
    }

    if (yy <= 0.5) {
        return(1.0 - 6.0 * yy^2 + 6.0 * yy^3)
    }

    2 * (1 - yy)^3
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Function for the Bartlett kernel.
### See Andrews, D, Econometrica, 1991, P821

kernelBartlett <- function (xx)
{
    yy <- abs(xx)
    (1 - yy) * as.numeric(yy <= 1.0)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### computes the \sigma^2_{mono} estimator, one could turn the
### is.montone flag off to get the \sigma^2_{pos} estimator.
### 
### For notation and details, see Geyer, Statistical Science, 1992,
### p473-483

MCObj <-
    function (statVals,
              verboseLevel = 0,
              lagCapProp   = 0.90,
              kern         = kernelBartlett,
              eps          = 0,
              isMonotone   = TRUE)
{
    stopifnot(is.function(kern))    
    stopifnot(lagCapProp >= 0.0)
    stopifnot(lagCapProp < 1.0)
    nn        <- length(statVals)
    lagCap    <- floor(nn * lagCapProp)
    llPrev    <- 2 * 1
    llNew     <- floor(min((nn - 1), 10 * sqrt(nn)))    

    ## making llNew even
    if (llNew %% 2 == 1) {
        llNew <- llNew + 1
    }
    
    count     <- 1
    sum       <- llNew
    flag      <- TRUE

    while (flag) {
        ptm <- proc.time( )
        if (verboseLevel >= 3) {
            cat('Note: [', count, '] In MCObj using lag =', llNew, ' ')
        }
        
        gamma <- acf(statVals, type = 'covariance', lag.max = llNew, plot = FALSE)
        gamma <- gamma$acf[ , 1, 1]
        rr    <- seq_len(llNew / 2)
        Gamma <- c(sapply(rr, FUN =
                          function (ii) { gamma[2 * ii - 1] + gamma[2 * ii] }))
        temp  <- rr[Gamma <= eps]

        time.used <- proc.time( ) - ptm
        if (verboseLevel >= 3) {       
            cat('[gamma:', round(gamma[llNew], 4),
                ', time:', round(time.used[1], 4), ']\n')       
        }
        
        if (length(temp) == 0) {
            count  <- count + 1
            llPrev <- llNew
            llNew  <- min(2 * llNew, lagCap)
            if (llNew %% 2 == 1) {
                llNew <- llNew + 1
            }
            sum <- sum + llNew
            next
        }

        flag  <- FALSE
        index <- min(temp) - 1
        rr    <- seq_len(index)
        stopifnot(all(Gamma[rr] > 0.0))
        if (isMonotone) {
            varEst <- -gamma[1] + 2 * sum(cummin(Gamma[rr]))
        }
        else {
            index  <- 2 * index
            varEst <- gamma[1] + (2 * sum(kern(seq.int(1, index - 1) / index) *
                                          gamma[seq.int(2, index)]))
        }
        if (varEst <= 0.0) {
            stop(paste('variance estimate in MCObj is non-positive:', varEst))
        }
    }
    
    if(verboseLevel >= 3) {
        cat('Note: In MCObj ', sum, ' many acf evaluations were used\n',
            'Note: In MCObj ', index, ' many lags were used\n', sep = '')
    }
    
    list(est       = mean(statVals),
         var       = varEst / nn,
         ESS       = nn / (varEst / gamma[1]),
         ESSFactor = varEst / gamma[1])
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ISObj <-
    function (statVals,
              lWeights,
              verboseLevel      = 0,
              kern              = kernelBartlett,
              lagCapProp        = 0.90,
              eps               = 0,
              useAdjustedWeight = FALSE,
              isMonotone        = FALSE)
{
    nn          <- length(statVals)
    weights     <- exp(lWeights)
    CC          <- mean(weights)
    vals        <- statVals * weights
    CVSqWeights <- var(weights) / CC^2
    list(est         = mean(vals) / CC,
         var         = MCObj(vals)$var / CC^2,
         ESS         = nn / (1 + CVSqWeights),
         CVSqWeights = CVSqWeights)
}

.get.postBurnIn <-
    function (EMCOut)
{
    nIters <- nrow(EMCOut$draws[ , , 1])
    burnIn <- min(floor(nIters / 4), 20000)
    seq.int((burnIn + 1), nIters)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Note the indirection through the use of getEMCOutFunc below in
### .findMaxTemper.do.1 and .placeTempers.do.4 is done to facilitate
### the package EMCC use the functionality developed in this
### package. Look at the use of EMC:::findMaxTemper and
### EMC:::placeTempers in EMCC/R/ladder.R

### Some helper functions for findMaxTemper

.check.acceptRatios <-
    function (acceptRatios, temper,
              acceptRatioTooLow  = 0.1,
              acceptRatioTooHigh = 0.7)
{
    if (acceptRatios[1] <= acceptRatioTooLow)
        cat('If possible, please make the "proposal density" for temperature = ',
            temper, '\n', 
            'less dispersed, the minimum observed acceptance rate (= ',
            acceptRatios[1], ') is below ', acceptRatioTooLow, '\n', sep = '')
    
    if (acceptRatios[1] >= acceptRatioTooHigh)
        cat('If possible, please make the "proposal density" for temperature = ',
            temper, '\n',
            'more dispersed, the minimum observed acceptance rate (= ',
            acceptRatios[1], ') is above ', acceptRatioTooHigh, '\n', sep = '')
}

.info.find.1 <-
    function (guideMe,
              acceptRatios,
              temper,
              detailedAcceptRatios,
              verboseLevel)
{
    if (!guideMe)
        return(invisible(0))
    
    .check.acceptRatios(acceptRatios, temper)    
}

.info.find.2 <-
    function (curr,
              prev,
              verboseLevel)
{
    if (verboseLevel >= 2) {
        cat('The computational steps for the provided statsFuncList:\n')
        cat('\n    MC estimates:', toString(round(curr$MCEsts, 4))) 
        cat('\n    IS estimates:', toString(round(curr$ISEsts, 4)))
        cat('\nMC var estimates:', toString(round(curr$MCVarEsts, 4))) 
        cat('\nIS var estimates:', toString(round(curr$ISVarEsts, 4))) 
        cat('\n           MCESS:', toString(round(curr$MCESS, 4))) 
        cat('\n           ISESS:', toString(round(curr$ISESS, 4)))
        cat('\n    DStat values:', toString(round(curr$DStats, 4)), '\n')
    }

    if (verboseLevel >= 3) {
        cat('\nprev fitness value summary:\n')
        print(summary(prev$fitness))
        cat('curr fitness value summary:\n')
        print(summary(curr$fitness))
        cat('prev statVals summary:\n')
        print(summary(prev$statVals))
        cat('curr statVals summary:\n')
        print(summary(curr$statVals))
        cat('curr logWeights summary:\n')
        print(summary(-prev$fitness * (1 / curr$temper - 1 / prev$temper)))
        cat('\n')
    }
}

.info.find.3 <-
    function (guideMe,
              curr,
              prev,
              cutoffESS,
              pos,
              verboseLevel)
{
    if (guideMe) {
        if (any(!is.na(curr$MCESS) &
                (curr$MCESS < cutoffESS)))
            warning('Please increase the sample size or improve the ',
                    '"proposal density";\n',
                    'some of the MCESSs (',  toString(round(curr$MCESS, 3)),
                    ') are too low\n',
                    call. = FALSE, immediate. = TRUE)
    
        if (any(!is.na(curr$ISESS) &
                (curr$ISESS < cutoffESS)))
            warning('Please increase the sample size or improve the ',
                    '"proposal density";\n',
                    'some of the ISESSs (', toString(round(curr$ISESS, 3)),
                    ') are too low\n',
                    call. = FALSE, immediate. = TRUE)
        
        .check.acceptRatios(curr$acceptRatios, curr$temper)
    }

    .info.find.2(curr         = curr,
                 prev         = prev,
                 verboseLevel = verboseLevel)
}

.get.infoListComp <-
    function (EMCOut,
              pbiSubset,
              EMCOutPos,
              statsFuncList,
              temper,
              guideMe,
              verboseLevel)
{
    detailedAcceptRatios <- EMCOut$detailedAcceptRatios
    acceptRatios         <- unlist(detailedAcceptRatios$MH[EMCOutPos, c(2, 4)])
    .info.find.1(guideMe              = guideMe,
                 acceptRatios         = acceptRatios,
                 temper               = temper,
                 detailedAcceptRatios = detailedAcceptRatios,
                 verboseLevel         = verboseLevel)

    nStats     <- length(statsFuncList)
    draws      <- EMCOut$draws[pbiSubset, , EMCOutPos]
    fitnessCol <- ncol(draws)

    ## throw away the fitness column use the rest of the draws matrix
    ## BUG REMOVED: in the previous version, we were not throwing away
    ## the fitness column
    
    statVals <- apply(draws[ , -fitnessCol], 1, function (row)
                      c(sapply(seq_along(statsFuncList), function (ss)
                               (statsFuncList[[ss]])(row))))
    
    if (length(statsFuncList) == 1)
        statVals <- matrix(statVals, ncol = 1)
    else
        statVals <- t(statVals)
    
    list(temper       = temper,
         draws        = draws,
         fitness      = draws[ , fitnessCol],
         statVals     = statVals,
         MCEsts       = numeric(nStats),
         MCVarEsts    = numeric(nStats),
         MCESS        = numeric(nStats),
         ISEsts       = numeric(nStats),
         ISVarEsts    = numeric(nStats),
         ISESS        = NA,
         DStats       = numeric(nStats),
         acceptRatios = acceptRatios)         
}

.update.infoListComp <-
    function (pos,
              curr,
              prev,
              statsFuncList,
              verboseLevel,
              eps = 0)
{
    if (verboseLevel >= 1) {
        cat('\n', rep('=', 50), '\n', sep = '')
        cat('Comparing temperatures [positions:', c(pos - 1, pos), ']:',
            toString(round(c(prev$temper, curr$temper), 4)), '\n')
    }

    logWeights <- -prev$fitness * (1 / curr$temper - 1 / prev$temper)
    for (ss in seq_along(curr$DStats)) {
        if (verboseLevel >= 2) {
            cat('Processing statistic:', ss, '\n')
        }
        
        nErrors <- 0

        ## the MC variance computation on the curr        
        tmpMC <- try(MCObj(curr$statVals[ , ss],
                           verboseLevel = verboseLevel,
                           eps          = eps),
                     silent = TRUE)
        if (inherits(tmpMC, 'try-error')) {
            nErrors <- nErrors + 1
            cat('\nThere was an error in computing the MC estimators for the ',
                ss, '-th statistic in the statsFuncList:\n', sep = '')
            print(statsFuncList[[ss]])
            cat('\nThe summary of the statistics values:\n')
            print(summary(curr$statVals[ , ss]))
            cat('\n')
            
            curr$MCEsts[ss]    <- NA
            curr$MCVarEsts[ss] <- NA
            curr$MCESS[ss]     <- NA
        }
        else { 
            curr$MCEsts[ss]    <- tmpMC$est
            curr$MCVarEsts[ss] <- tmpMC$var
            curr$MCESS[ss]     <- tmpMC$ESS
        }
        
        ## the IS variance computation on the prev        
        tmpIS <- try(ISObj(prev$statVals[ , ss],
                           lWeights     = logWeights,
                           verboseLevel = verboseLevel,
                           eps          = eps),
                     silent = TRUE)
        if (inherits(tmpIS, 'try-error')) {
            nErrors <- nErrors + 1
            cat('\nThere was an error in computing the IS estimators for the ',
                ss, '-th statistic in the statsFuncList:\n', sep = '')
            print(statsFuncList[[ss]])
            cat('\nThe summary of the statistics values:\n')
            print(summary(prev$statVals[ , ss]))
            cat('\n')
            
            curr$ISEsts[ss]    <- NA
            curr$ISVarEsts[ss] <- NA
            curr$ISESS         <- NA
        }
        else {         
            curr$ISEsts[ss]    <- tmpIS$est
            curr$ISVarEsts[ss] <- tmpIS$var
            curr$ISESS         <- tmpIS$ESS
        }

        if (nErrors > 0)
            curr$DStats[ss] <- NA
        else
            curr$DStats[ss]    <- ((tmpMC$est - tmpIS$est) /
                                   sqrt(tmpMC$var + tmpIS$var))        
    }
    
    curr
}

.findMaxTemper.do.2 <-
    function (startLadder,
              infoList,
              nStats,
              cutoffDStats,
              pbiSubset,
              sampDim,
              levelsSaveSampFor,
              saveFitness,
              verboseLevel)
{
    nLevels <- length(startLadder)
    ret     <- list(temperLadder      = startLadder,
                    MCEsts            = matrix(0, nLevels - 1, nStats),
                    MCVarEsts         = matrix(0, nLevels - 1, nStats),
                    MCESS             = matrix(0, nLevels - 1, nStats),
                    ISEsts            = matrix(0, nLevels - 1, nStats),
                    ISVarEsts         = matrix(0, nLevels - 1, nStats),
                    ISESS             = numeric(nLevels - 1),
                    DStats            = matrix(0, nLevels - 1, nStats),
                    cutoffDStats      = cutoffDStats,
                    nIters            = length(pbiSubset),
                    levelsSaveSampFor = levelsSaveSampFor,
                    draws             = NA)

    for (ii in seq_len(nLevels - 1)) {
        ret$MCEsts[ii, ]    <- infoList[[ii + 1]]$MCEsts
        ret$MCVarEsts[ii, ] <- infoList[[ii + 1]]$MCVarEsts
        ret$MCESS[ii, ]     <- infoList[[ii + 1]]$MCESS
        ret$ISEsts[ii, ]    <- infoList[[ii + 1]]$ISEsts
        ret$ISVarEsts[ii, ] <- infoList[[ii + 1]]$ISVarEsts
        ret$ISESS[ii]       <- infoList[[ii + 1]]$ISESS
        ret$DStats[ii, ]    <- infoList[[ii + 1]]$DStats
    }

    nInfoList        <- length(infoList)    
    ret$temperLadder <- numeric(nInfoList)
    for (ii in seq_len(nInfoList)) {
        ret$temperLadder[ii] <- infoList[[ii]]$temper
    }

    if ((levelsSaveSampForLen <- length(levelsSaveSampFor)) > 0) {
        dd        <- ifelse(saveFitness, sampDim + 1, sampDim)
        tmp       <- seq_len(dd)
        ret$draws <- array(dim = c(ret$nIters, dd, levelsSaveSampForLen))
        for (ii in seq_len(levelsSaveSampForLen)) {
            jj <- levelsSaveSampFor[ii]
            if (is.null(infoList[[jj]])) {
                warning('Could not save samples for level: ', jj,
                        'since sampling was not done for it',
                        call. = FALSE, immediate. = TRUE)
                next
            }
            ret$draws[ , , ii] <- infoList[[jj]]$draws[ , tmp]
        }
    }        
    ret
}

.findMaxTemper.do.1 <-
    function (getEMCOutFunc,
              statsFuncList,
              temperLadder,
              cutoffDStats,
              cutoffESS,
              guideMe,
              levelsSaveSampFor,
              saveFitness,
              doFullAnal,
              verboseLevel,
              ...)
{
    if (!doFullAnal) startLadder <- temperLadder[c(1, 2)]
    else             startLadder <- temperLadder

    nStats  <- length(statsFuncList)
    nLevels <- length(temperLadder)    
    
    ## whether a plausible maximum was foundForAll
    foundForAll <- FALSE
    toFindFor   <- seq_len(nStats)
    foundSoFar  <- c( )
    
    ## whether scanned all tempers
    scannedAll <- FALSE
    pos        <- 1

    while ((!foundForAll) && (!scannedAll)) {
        if (pos == 1) {
            EMCOut    <- getEMCOutFunc(startLadder, ...)
            sampDim   <- ncol(EMCOut$draws[ , , 1]) - 1
            pbiSubset <- .get.postBurnIn(EMCOut)
            EMCOutPos <- pos
            infoList  <- list( )
        }
       
        infoList[[pos]] <- .get.infoListComp(EMCOut,
                                             pbiSubset     = pbiSubset,
                                             EMCOutPos     = EMCOutPos,
                                             statsFuncList = statsFuncList,
                                             temper        = startLadder[pos],
                                             guideMe       = guideMe,
                                             verboseLevel  = verboseLevel)
        if (pos == 1) {
            pos       <- pos + 1
            EMCOutPos <- pos
            next
        }

        curr            <- infoList[[pos]]
        prev            <- infoList[[pos - 1]]
        curr            <- .update.infoListComp(pos           = pos,
                                                curr          = curr,
                                                prev          = prev,
                                                statsFuncList = statsFuncList,
                                                verboseLevel  = verboseLevel)
        infoList[[pos]] <- curr
        
        .info.find.3(guideMe      = guideMe,
                     curr         = curr,
                     prev         = prev,
                     cutoffESS    = cutoffESS,
                     pos          = pos,
                     verboseLevel = verboseLevel)
        
        foundFor <- toFindFor[!is.na(curr$DStats) &
                              (abs(curr$DStats) > cutoffDStats)]
        if (length(foundFor) > 0) {
            if (verboseLevel >= 1) {
                cat('Found a plausible max temper: ', prev$temper, ', for ',
                    'statistic(s): ', toString(foundFor), '\n', sep = '')
            }

            if (!doFullAnal) {
                foundSoFar <- union(foundSoFar, foundFor)

                if (verboseLevel >= 1) {
                    cat('Found plausible temperatures so far, for statistic(s):',
                        toString(foundSoFar), '\n\n')
                }

                ## found plausible temperatures for all the
                ## statistics, so stop                
                if ((length(setdiff(toFindFor, foundSoFar)) == 0)) {
                    levelsSaveSampFor <- seq_len(pos)
                    foundForAll       <- TRUE
                    next
                }
            }
        }

        ## did not find max temperature, the ladder was exhausted, so
        ## stop
        if (pos == nLevels) {
            scannedAll <- TRUE
            next
        }
        
        pos       <- pos + 1
        EMCOutPos <- pos
        if (doFullAnal) {
            next
        }
        
        ## !doFullAnal, then need to update the infoList by
        ## introducing another ladder in the startLadder
        startLadder <- temperLadder[seq_len(pos)]
        EMCOut      <- getEMCOutFunc(startLadder[pos], ...)
        EMCOutPos   <- 1
    }

    .findMaxTemper.do.2(startLadder       = startLadder,
                        infoList          = infoList,
                        nStats            = nStats,
                        cutoffDStats      = cutoffDStats,
                        pbiSubset         = pbiSubset,
                        sampDim           = sampDim,
                        levelsSaveSampFor = levelsSaveSampFor,
                        saveFitness       = saveFitness,
                        verboseLevel      = verboseLevel)
}

.findMaxTemper.do <-
    function (nIters,
              statsFuncList,
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc,
              temperLadder,
              cutoffDStats,
              cutoffESS,
              guideMe,
              levelsSaveSampFor,
              saveFitness,
              doFullAnal,
              verboseLevel,
              getEMCOutFunc = NULL,
              ...)
{
    if (is.null(getEMCOutFunc)) {
        getEMCOutFunc <-
            function (ladder, ...)
            {
                EMCOut <- 
                    TOEMCMain(nIters            = nIters,
                              temperLadder      = ladder,
                              startingVals      = startingVals,
                              logTarDensFunc    = logTarDensFunc,
                              MHPropNewFunc     = MHPropNewFunc, 
                              logMHPropDensFunc = logMHPropDensFunc,
                              moveProbsList     = list(MH = 1.0),
                              moveNTimesList    = list(MH = 1),
                              levelsSaveSampFor = seq_along(ladder),
                              saveFitness       = TRUE,
                              verboseLevel      = verboseLevel,
                              ...)
                
                if (verboseLevel >= 1) {
                    cat('\nThe detailed MH acceptance ratios:\n')
                    print(EMCOut$detailedAcceptRatios$MH)
                }
                
                EMCOut
            }
    }

    .findMaxTemper.do.1(getEMCOutFunc     = getEMCOutFunc,
                        statsFuncList     = statsFuncList,
                        temperLadder      = temperLadder,
                        cutoffDStats      = cutoffDStats,
                        cutoffESS         = cutoffESS,
                        guideMe           = guideMe,
                        levelsSaveSampFor = levelsSaveSampFor,
                        saveFitness       = saveFitness,
                        doFullAnal        = doFullAnal,
                        verboseLevel      = verboseLevel,
                        ...)
}

findMaxTemper <-
    function (nIters,
              statsFuncList,
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc = NULL, 
              temperLadder      = NULL,
              temperLimits      = NULL,
              ladderLen         = 10,
              scheme            = 'exponential',
              schemeParam       = 0.5,
              cutoffDStats      = 1.96,
              cutoffESS         = 50,
              guideMe           = TRUE,
              levelsSaveSampFor = NULL,
              saveFitness       = FALSE,
              doFullAnal        = TRUE,
              verboseLevel      = 0,
              ...)   
{
    ptm <- proc.time( )
    ## Error checks
    statsFuncList <- .check.statsFuncList(statsFuncList)
    temperLadder  <- .get.temperLadder(temperLadder = temperLadder,
                                       temperLimits = temperLimits,
                                       ladderLen    = ladderLen,
                                       scheme       = scheme,
                                       schemeParam  = schemeParam)
    nLevels       <- length(temperLadder)
    cutoffDStats  <- .check.numericWithLLim(cutoffDStats, 0, retFunc = as.numeric)
    cutoffESS     <- .check.numericWithLLim(cutoffESS, 50)
    guideMe       <- .check.logical(guideMe)    

    levelsSaveSampFor <- .check.levelsSaveSamplesFor(levelsSaveSampFor, nLevels)
    saveFitness       <- .check.logical(saveFitness)
    doFullAnal        <- .check.logical(doFullAnal)
    verboseLevel      <- .check.numericWithLLim(verboseLevel, NA)

    if (verboseLevel >= 1) cat('\nBEGIN: find maximum temperature\n')    
    ret <- .findMaxTemper.do(nIters            = nIters,
                             statsFuncList     = statsFuncList,
                             startingVals      = startingVals,
                             logTarDensFunc    = logTarDensFunc,
                             MHPropNewFunc     = MHPropNewFunc,
                             logMHPropDensFunc = logMHPropDensFunc,
                             temperLadder      = temperLadder,
                             cutoffDStats      = cutoffDStats,
                             cutoffESS         = cutoffESS,
                             guideMe           = guideMe,
                             levelsSaveSampFor = levelsSaveSampFor,
                             saveFitness       = saveFitness,
                             doFullAnal        = doFullAnal,
                             verboseLevel      = verboseLevel,
                             ...)
    if (verboseLevel >= 1) cat('E N D: find maximum temperature\n')    
    if (!any(is.na(ptm))) ret$time <- proc.time( ) - ptm
    ret$startingVals <- startingVals
    class(ret)       <- 'EMCMaxTemper'
    ret
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.EMCMaxTemper <-
    function (x, ...)
{
    dd       <- dim(x$DStats)
    rr1      <- seq_len(dd[1])
    rr2      <- seq_len(dd[2])
    rr2names <- paste('statistic', rr2, sep = '')

    cat('\nThe temperature and the corresponding D-statistics for',
        'various statistics:\n\n')
    tmp           <- as.matrix(cbind(x$temperLadder,
                                     rbind(rep(NA, dd[2]), x$DStats)))
    dimnames(tmp) <-
        list(paste('level', seq_along(x$temperLadder), sep = ''),
             c('temperature', rr2names))
    print(tmp, ...)

    cat('\nSuggested maximum temperatures by various statistics ',
        'with cutoff = ', x$cutoffDStats, ':\n\n', sep = '')
    tmp    <- data.frame(level       = rr2,
                         temperature = NA)
    for (ii in rr2) {
        filter <- (!is.na(x$DStats[ , ii]) &
                   abs(x$DStats[ , ii]) >= x$cutoffDStats)
        if (any(filter)) {
            jj         <- min(rr1[filter])
            tmp[ii, 1] <- jj
            tmp[ii, 2] <- as.character(round(x$temperLadder[jj], 4))
        }
        else {
            jj         <- dd[1] + 1
            tmp[ii, 1] <- jj
            tmp[ii, 2] <- paste(round(x$temperLadder[jj], 4),
                                '[all below cut-off]')
        }
    }
    rownames(tmp) <- rr2names
    print(tmp)
    cat('\n')
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Some helper functions for placeTempers

ARHatIS <-
    function (fitness1,
              fitness2,
              tau1,
              tau2,
              temper1,
              temper2)
{
    invTau1    <- 1.0 / tau1
    invTau2    <- 1.0 / tau2
    invTemper1 <- 1.0 / temper1
    invTemper2 <- 1.0 / temper2

    HVals          <- exp(pmin(0, (fitness1 - fitness2) *
                               (invTemper1 - invTemper2)))
    lWeights       <- -(fitness1 * (invTemper1 - invTau1) +
                        fitness2 * (invTemper2 - invTau2))
    adjWeights     <- exp(lWeights - max(lWeights))
    meanAdjWeights <- mean(adjWeights)
        
    list(est         = mean(HVals * adjWeights) / meanAdjWeights,
         CVSqWeights = var(adjWeights) / (meanAdjWeights^2))
}

.info.place.1 <-
    function (guideMe,
              detailedAcceptRatios,
              temperLadder,
              verboseLevel)
{
    if (!guideMe)
        return(invisible(0))

    acceptRatiosMH <- detailedAcceptRatios$MH
    for (ii in seq_along(temperLadder)) 
        .check.acceptRatios(unlist(acceptRatiosMH[ii, c(2, 4)]), temperLadder[ii])
    
    if (verboseLevel >= 1) {
        cat('The detailed MH accept ratios:\n')
        print(acceptRatiosMH)
        cat('The detailed RE accept ratios:\n')
        print(detailedAcceptRatios$BetweenAllLevels)
    }    
}

.info.place.2 <-
    function (startLadder,
              currInterval,
              verboseLevel)
{
    if (verboseLevel >= 2) {
        startLadderLen <- length(startLadder)
        cat(rep('o', 80), '\n', 
            'BEGIN: Trying with ', startLadderLen, ' temperatures...\n',
            'This temperature ladder [length = ', startLadderLen, ']:\n',
            sep = '')
        print(startLadder)
    }

    if (verboseLevel >= 1) {
        cat(rep('=', 50), '\n', sep = '')        
        cat('Trying to place temperature for interval:',
            toString(round(currInterval, 4)), '\n\n')
    }
}        

.info.place.3 <-
    function (nShrink,
              nExpand,
              propInterval,
              acceptRatioEst,
              verboseLevel)
{
    if (verboseLevel >= 1) {
        tmp           <- data.frame(nShrink          = nShrink,
                                    nExpand          = nExpand,
                                    acceptedInterval = toString(round(propInterval, 4)),
                                    acceptRatioEst   = acceptRatioEst)
        rownames(tmp) <- 'Summary of trials:'
        print(tmp)
    }
}

.info.place.4 <-
    function (startLadder,
              endString,
              CVSqWeights,
              ladderSmall,
              badAceeptRatioLimits,
              guideMe,
              verboseLevel)
{
    if (verboseLevel >= 2) {
        cat(rep('o', 80), '\n',
            'The minimum number of temperatures needed:', length(startLadder), '\n',
            'NOTE: ', endString, '\n',
            'NOTE: The summary of CVSq of weights is:\n', sep = '')
        print(summary(CVSqWeights))
    }

    if (guideMe) {
        if (ladderSmall)
            cat(endString, '\n')
    }
    
    if (badAceeptRatioLimits)
        cat('[NOTE:', endString, ']\n')        
}

.placeTempers.do.4 <-
    function (nIters,
              acceptRatioLimits,
              ladderLenMax,
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc,
              temperLadder,
              temperLimits,
              ladderLen,
              guideMe,
              levelsSaveSampFor,
              saveFitness,
              verboseLevel,
              getEMCOutFunc = NULL,
              ...)
{
    if (!is.null(getEMCOutFunc))
        return(getEMCOutFunc(temperLadder, ...))
    
    nLevels <- length(temperLadder)
    EMCOut  <- TOEMCMain(nIters            = nIters,
                         temperLadder      = temperLadder,
                         startingVals      = startingVals,
                         logTarDensFunc    = logTarDensFunc,
                         MHPropNewFunc     = MHPropNewFunc, 
                         logMHPropDensFunc = logMHPropDensFunc,
                         moveProbsList     = list(MH = 0.3, RE = 0.7),
                         moveNTimesList    = list(MH = 1, RE = nLevels),
                         levelsSaveSampFor = seq_len(nLevels),
                         saveFitness       = TRUE,
                         verboseLevel      = verboseLevel,
                         ...)
    EMCOut
}    

.placeTempers.do.3 <-
    function (EMCOut,
              startLadder,
              ARHats,
              CVSqWeights,
              temperLadder,
              temperLimits,
              acceptRatioLimits,
              pbiSubset,
              sampDim,
              levelsSaveSampFor,
              saveFitness,
              verboseLevel)
{
    ret <- list(finalLadder      = startLadder,
                temperLadder      = temperLadder,
                acceptRatiosEst   = ARHats,
                CVSqWeights       = CVSqWeights,
                temperLimits      = sort(temperLimits, decreasing = TRUE),
                acceptRatioLimits = acceptRatioLimits,
                nIters            = length(pbiSubset),
                levelsSaveSampFor = levelsSaveSampFor,
                draws             = NA)

    if ((levelsSaveSampForLen <- length(levelsSaveSampFor)) > 0) {
        dd        <- ifelse(saveFitness, sampDim + 1, sampDim)
        tmp       <- seq_len(dd)
        ret$draws <- array(dim = c(ret$nIters, dd, levelsSaveSampForLen))
        for (ii in seq_len(levelsSaveSampForLen)) {
            jj                 <- levelsSaveSampFor[ii]
            ret$draws[ , , ii] <- EMCOut$draws[pbiSubset, tmp, jj]
        }
    }
    ret    
}

.placeTempers.do.2 <-
    function (EMCOut,
              temperLadder,
              propInterval,
              pbiSubset,
              sampDim,
              verboseLevel)
{
    temper1     <- propInterval[1]
    temper2     <- propInterval[2]
    tmp         <- seq_along(temperLadder)
    jj1         <- max(tmp[temperLadder >= temper1])
    jj2         <- max(tmp[temperLadder >= temper2])
    tau1        <- temperLadder[jj1]
    tau2        <- temperLadder[jj2]
    fitness1    <- EMCOut$draws[pbiSubset, sampDim + 1, jj1]
    fitness2    <- EMCOut$draws[pbiSubset, sampDim + 1, jj2]
    ARHatISList <- ARHatIS(fitness1 = fitness1,
                           fitness2 = fitness2,
                           tau1     = tau1,
                           tau2     = tau2,
                           temper1  = temper1,
                           temper2  = temper2)

    if (verboseLevel >= 2) {
        cat('\nimportance sampling estimate of acceptRatio:',
            round(ARHatISList$est, 4))
        cat('\n                      for proposed interval:',
            toString(round(c(temper1, temper2), 4)))
        cat('\n                         using temperatures:',
            toString(round(c(tau1, tau2), 4)))
        cat('\n                       with CVSq of weights:',
            round(ARHatISList$CVSqWeights, 4), '\n')
    }

    if (abs(tau1 - tau2) <= 0) {
        msg <- paste('The provided temperature ladder is too small, some more',
                     'intermediate temperatures are needed in the interval: ',
                     toString(round(c(temper1, temper2), 4)))
        fatal(msg)
    }
    ARHatISList
}

.placeTempers.do.1 <-
    function (EMCOut,
              acceptRatioLimits,
              ladderLenMax,
              temperLadder,
              temperLimits,
              ladderLen,
              guideMe,
              levelsSaveSampFor,
              saveFitness,
              verboseLevel)
{
    pbiSubset          <- .get.postBurnIn(EMCOut)
    sampDim            <- ncol(EMCOut$draws[ , , 1]) - 1

    ulimTL             <- temperLadder[1]
    llimTL             <- temperLadder[length(temperLadder)]
    startLadder        <- c(ulimTL, llimTL)

    acceptRatioLimits  <- sort(acceptRatioLimits, decreasing = FALSE)
    ulimAR             <- acceptRatioLimits[2]
    llimAR             <- acceptRatioLimits[1]

    ARComputed         <- c(FALSE)
    ARHatsWithinLimits <- c(FALSE)
    ARHats             <- c(0)
    CVSqWeights        <- c(0)
    ladderSmall        <- FALSE
    found              <- FALSE
    ii                 <- 1

    while (!found) {
        startLadderLen <- length(startLadder)
        currInterval   <- startLadder[c((startLadderLen - 1), startLadderLen)]
        propInterval   <- currInterval

        .info.place.2(startLadder  = startLadder,
                      currInterval = currInterval,
                      verboseLevel = verboseLevel)

        insertedLoopCount    <- 0
        insertedLoopCountMax <- 100
        nShrink              <- 0
        nExpand              <- 0
        badAceeptRatioLimits <- FALSE        
        inserted             <- FALSE
        changed              <- TRUE
        while ((!inserted) && changed) {
            insertedLoopCount <- insertedLoopCount + 1
            ARHatISList      <- .placeTempers.do.2(EMCOut       = EMCOut,
                                                   temperLadder = temperLadder,
                                                   propInterval = propInterval,
                                                   pbiSubset    = pbiSubset,
                                                   sampDim      = sampDim,
                                                   verboseLevel = verboseLevel)
            ARHats[ii]        <- ARHatISList$est

            ## ARHat too low: shrinking interval
            if (ARHats[ii] < llimAR) {
                nShrink         <- nShrink + 1
                ll              <- propInterval[2]
                uu              <- propInterval[1]                    
                propInterval[2] <- runif(1, ll, mean(c(ll, uu)))
                if (verboseLevel >= 2) {
                    cat('Shrinked interval [with acceptance ratio: ',
                        ARHats[ii], ']:\n', sep = '')
                    print(propInterval)
                }
                next
            }

            ## ARHat too high: expanding interval
            if (ARHats[ii] > ulimAR) {
                nExpand         <- nExpand + 1
                ll              <- currInterval[2]
                uu              <- propInterval[2]
                propInterval[2] <- runif(1, ll, mean(c(uu, ll)))
                if (abs((propInterval[2] - currInterval[2]) <= 0.0)) {
                    inserted <- TRUE
                }

                if (verboseLevel >= 2) {
                    cat('Expanded interval [with acceptRatio: ',
                        ARHats[ii], ']:\n', sep = '')
                    print(propInterval)
                }
                next
            }
            
            ## ARHat moderate: (llimAR <= ARHats[ii]) && (ARHats[ii] <= ulimAR)
            ARHatsWithinLimits[ii] <- TRUE
            CVSqWeights[ii]        <- ARHatISList$CVSqWeights
            CVSqWeights            <- c(CVSqWeights, 0)
            
            ## expand the ladder
            if (abs(propInterval[2] - currInterval[2]) > 0.0) {
                ARHats             <- c(ARHats, 0)
                ARHatsWithinLimits <- c(ARHatsWithinLimits, FALSE)
                startLadder        <- c(startLadder[seq.int(1, ii)],
                                        propInterval[2],
                                        startLadder[seq.int(ii + 1, startLadderLen)])
                ii                 <- ii + 1
                inserted           <- TRUE
                thisII             <- ii - 1
            }
            else {
                changed            <- FALSE
                thisII             <- ii
            }

            if (insertedLoopCount > insertedLoopCountMax) {
                badAceeptRatioLimits <- TRUE
                found                <- TRUE
                endString            <-
                    paste('Please change the acceptRatioLimits and re-run',
                          'this function; the temperature ladder provided',
                          'as the output is the log-scale placement, it may not',
                          'be the "best" or the "shortest"')
                startLadder          <- buildLadder(limits = temperLimits,
                                                    length = ladderLen,
                                                    scheme = 'log')
                ARHats               <- NULL
                CVSqWeights          <- NULL
            }

            .info.place.3(nShrink        = nShrink,
                          nExpand        = nExpand,
                          propInterval   = propInterval,
                          acceptRatioEst = ARHats[thisII],
                          verboseLevel   = verboseLevel)
        }

        if (abs((propInterval[2] - currInterval[2]) <= 0.0)) {
            found     <- TRUE
            endString <- 'The acceptRatioLimits criterion was achieved'
        }        

        if (ladderLenMax <= startLadderLen) {            
            found       <- TRUE
            ladderSmall <- TRUE
            endString   <- paste('The ladderLenMax criterion was achieved,',
                                 'you may want to increase it')
        }

        if (verboseLevel >= 2) {
            cat('The acceptance ratios so far:\n')
            print(ARHats)
            cat('The acceptance ratio within acceptRatioLimits indicator:\n')
            print(ARHatsWithinLimits)        
            cat('END: Trying with', startLadderLen, 'temperatures...\n')
        }                        
    }

    .info.place.4(startLadder          = startLadder,
                  endString            = endString,
                  CVSqWeights          = CVSqWeights,
                  ladderSmall          = ladderSmall,
                  badAceeptRatioLimits = badAceeptRatioLimits,
                  guideMe              = guideMe,
                  verboseLevel         = verboseLevel)

    .placeTempers.do.3(EMCOut            = EMCOut,
                       startLadder       = startLadder,
                       ARHats            = ARHats,
                       CVSqWeights       = CVSqWeights,
                       temperLadder      = temperLadder,
                       temperLimits      = temperLimits,
                       acceptRatioLimits = acceptRatioLimits,
                       pbiSubset         = pbiSubset,
                       sampDim           = sampDim,
                       levelsSaveSampFor = levelsSaveSampFor,
                       saveFitness       = saveFitness,
                       verboseLevel      = verboseLevel)
}

.placeTempers.do <-
    function (nIters,
              acceptRatioLimits,
              ladderLenMax,
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc,
              temperLadder,
              temperLimits,
              ladderLen,
              guideMe,
              levelsSaveSampFor,
              saveFitness,
              verboseLevel,
              ...)
{
    EMCOut <- .placeTempers.do.4(nIters            = nIters,
                                 acceptRatioLimits = acceptRatioLimits,
                                 ladderLenMax      = ladderLenMax,
                                 startingVals      = startingVals,
                                 logTarDensFunc    = logTarDensFunc,
                                 MHPropNewFunc     = MHPropNewFunc,
                                 logMHPropDensFunc = logMHPropDensFunc,
                                 temperLadder      = temperLadder,
                                 temperLimits      = temperLimits,
                                 ladderLen         = ladderLen,
                                 guideMe           = guideMe,
                                 levelsSaveSampFor = levelsSaveSampFor,
                                 saveFitness       = saveFitness,
                                 verboseLevel      = verboseLevel,
                                 ...)

    .info.place.1(guideMe              = guideMe,
                  detailedAcceptRatios = EMCOut$detailedAcceptRatios,
                  temperLadder         = temperLadder,
                  verboseLevel         = verboseLevel)

    .placeTempers.do.1(EMCOut            = EMCOut,
                       acceptRatioLimits = acceptRatioLimits,
                       ladderLenMax      = ladderLenMax,
                       temperLadder      = temperLadder,
                       temperLimits      = temperLimits,
                       ladderLen         = ladderLen,
                       guideMe           = guideMe,
                       levelsSaveSampFor = levelsSaveSampFor,
                       saveFitness       = saveFitness,
                       verboseLevel      = verboseLevel)
}

placeTempers <-
    function (nIters,
              acceptRatioLimits,
              ladderLenMax,              
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc = NULL,              
              temperLadder      = NULL,
              temperLimits      = NULL,
              ladderLen         = 15,
              scheme            = 'exponential',
              schemeParam       = 1.5,
              guideMe           = TRUE,
              levelsSaveSampFor = NULL,
              saveFitness       = FALSE,
              verboseLevel      = 0,
              ...)                 
{
    ptm <- proc.time( )
    ## Error checks
    acceptRatioLimits <- .check.acceptRatioLimits(acceptRatioLimits)
    temperLadder      <- .get.temperLadder(temperLadder = temperLadder,
                                           temperLimits = temperLimits,
                                           ladderLen    = ladderLen,
                                           scheme       = scheme,
                                           schemeParam  = schemeParam)
    ladderLenMax      <- .check.numericWithLLim(ladderLenMax, 3)
    guideMe           <- .check.logical(guideMe)
    nLevels           <- length(temperLadder)

    levelsSaveSampFor <- .check.levelsSaveSamplesFor(levelsSaveSampFor, nLevels)
    saveFitness       <- .check.logical(saveFitness)
    verboseLevel      <- .check.numericWithLLim(verboseLevel, NA)

    if (verboseLevel >= 1) cat('\nBEGIN: place temperatures\n')    
    ret <- .placeTempers.do(nIters            = nIters,
                            acceptRatioLimits = acceptRatioLimits,
                            ladderLenMax      = ladderLenMax,
                            startingVals      = startingVals,
                            logTarDensFunc    = logTarDensFunc,
                            MHPropNewFunc     = MHPropNewFunc,
                            logMHPropDensFunc = logMHPropDensFunc,
                            temperLadder      = temperLadder,
                            temperLimits      = temperLimits,
                            ladderLen         = ladderLen,
                            guideMe           = guideMe,
                            levelsSaveSampFor = levelsSaveSampFor,
                            saveFitness       = saveFitness,
                            verboseLevel      = verboseLevel,
                            ...)
    if (verboseLevel >= 1) cat('E N D: place temperature\n')
    if (!any(is.na(ptm))) ret$time <- proc.time( ) - ptm
    ret$startingVals <- startingVals
    class(ret) <- 'EMCPlaceTempers'
    ret
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.EMCPlaceTempers <-
    function (x, ...)
{
    cat('\nThe suggested final temperature ladder for acceptance rate range [',
        toString(x$acceptRatioLimits), ']:\n', sep = '')
    print(x$finalLadder, ...)
    cat('\n')
}    
