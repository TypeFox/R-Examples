
###  $Id: doEMC.R,v 1.35 2008/07/06 03:09:12 goswami Exp $
###  
###  File:    doEMC.R
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
### The following is the main 'workhorse' function

TOEMCMain <-
    function (nIters,              
              temperLadder,
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc        = NULL, 
              MHRWMPropDisp            = NULL,
              MHBlocks                 = NULL,
              MHBlockNTimes            = NULL,
              isMHRWMVec               = NULL,
              moveProbsList            = NULL,
              moveNTimesList           = NULL,
              SCRWMNTimes              = NULL,
              SCRWMPropSD              = NULL,
              moveSelectionCodesList   = NULL,
              moveSelectionTempersList = NULL,
              levelsSaveSampFor        = NULL,
              nThin                    = 1,
              saveFitness              = FALSE,
              saveAcceptRatiosList     = FALSE, 
              timeInSecs               = -1,
              verboseLevel             = 0,
              ...)   
{
    ptm <- proc.time( )
    ## BEGIN: Error checks
    nIters       <- .check.numericWithLLim(nIters, 0)
    temperLadder <- .check.temperLadder(temperLadder)
    nLevels      <- as.integer(length(temperLadder))
    startingVals <- .check.startingVals(startingVals, nLevels)
    sampDim      <- as.integer(ncol(startingVals))
    
    logTarDensFunc    <- .check.logTarDensFunc(logTarDensFunc)
    MHPropNewFunc     <- .check.MHPropNewFunc(MHPropNewFunc)
    logMHPropDensFunc <- .check.logMHPropDensFunc(logMHPropDensFunc)
    oneIterMHFunc     <- .oneIterMHFuncGen(lTDF   = logTarDensFunc,
                                           MHPNF  = MHPropNewFunc,
                                           lMHPDF = logMHPropDensFunc)
    
    MHRWMPropDispArr     <- .check.MHRWMPropDisp(MHRWMPropDisp,
                                                 sampDim       = sampDim,
                                                 nLevels       = nLevels,
                                                 temperLadder  = temperLadder,
                                                 MHPropNewFunc = MHPropNewFunc)

    MHBlocks         <- .check.MHBlocks(MHBlocks, sampDim)
    MHNBlocks        <- as.integer(length(MHBlocks))
    MHBlockNTimes    <- .check.MHBlockNTimes(MHBlockNTimes, MHNBlocks)
    isMHRWMVec       <- .check.isMHRWMVec(isMHRWMVec, MHNBlocks)

    moveProbsList  <- .check.moveProbsList(moveProbsList, nLevels)
    moveNTimesList <- .check.moveNTimesList(moveNTimesList,
                                            moveProbsList = moveProbsList,
                                            nLevels       = nLevels)
    
    SCRWMPropSD <- .check.SCRWMPropSD(SCRWMPropSD, moveNTimesList)
    SCRWMNTimes <- .check.SCRWMNTimes(SCRWMNTimes, moveNTimesList)

    moveSelectionCodesList   <- .check.mSCList(moveSelectionCodesList,
                                               moveNTimesList = moveNTimesList)
    moveSelectionTempersList <- .check.mSTList(moveSelectionTempersList,
                                               moveNTimesList = moveNTimesList,
                                               temperColdest  = temperLadder[nLevels])

    levelsSaveSampFor <- .check.levelsSaveSamplesFor(levelsSaveSampFor, nLevels)
    nThin             <- .check.numericWithLLim(nThin, 1)
    saveFitness       <- .check.logical(saveFitness)    
    timeInSecs        <- .check.timeInSecs(timeInSecs)      
    verboseLevel      <- .check.numericWithLLim(verboseLevel, NA)  
    procTimeFunc      <- as.function(proc.time)             
    procTimeFuncEnv   <- new.env( )
    doCallFunc        <- as.function(doCall)             
    doCallFuncEnv     <- new.env( )
    dotsList          <- list(...)
    nSave             <- nIters %/% nThin
    if ((nIters %% nThin) != 0) {
        ## Increasing the nSave to make sure we have enough space to
        ## save all the iters
        nSave <- nSave + 1
        warning('nIters: ', nIters, ' is not divisible by nThin: ', nThin)
    }
    nSave    <- as.integer(nSave)
    argsList <- collectVarnames(ls( ))
    ## E N D: Error checks

    if (argsList$verboseLevel >= 3) {
        cat('The processed arguments:\n')
        print(argsList)
    }
    if (argsList$verboseLevel >= 1) cat('\nBEGIN: EMC\n')        
    res <- .Call('TOEMCMainC', argsList, PACKAGE = 'EMC')
    if (argsList$verboseLevel >= 1) cat('E N D: EMC\n')

    res$nIters               <- nIters
    res$nThin                <- nThin
    res$nSave                <- nSave
    res$levelsSaveSampFor    <- levelsSaveSampFor
    res$temperLadder         <- temperLadder
    res$startingVals         <- startingVals
    res$moveProbsList        <- moveProbsList
    res$moveNTimesList       <- moveNTimesList
    res$detailedAcceptRatios <- .get.detailedAcceptRatios(res$acceptRatiosList)

    if (!saveAcceptRatiosList)
        res$acceptRatiosList <- NULL

    if (!any(is.na(ptm))) res$time <- proc.time( ) - ptm
    class(res) <- 'EMC'
    res
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomWalkMetropolis <-
    function (nIters,              
              startingVal,
              logTarDensFunc,
              propNewFunc,
              MHBlocks      = NULL,
              MHBlockNTimes = NULL,
              nThin         = 1,
              saveFitness   = FALSE, 
              verboseLevel  = 0,
              ...)     
{
    if (is.null(propNewFunc)) {
        MHPropNewFunc <- NULL
    }
    else {
        propNewFunc <- .check.func.do(propNewFunc, 
                                      argsReq   = c('block', 'currentDraw', '...'),
                                      retObjMsg = 'a numeric vector')
        MHPropNewFunc <- 
            function (temperature, block, currentDraw, ...)
                propNewFunc(block, currentDraw, ...)     
    }

    res <- TOEMCMain(nIters            = nIters,              
                     temperLadder      = c(1),
                     startingVals      = startingVal,
                     logTarDensFunc    = logTarDensFunc,
                     logMHPropDensFunc = NULL, 
                     MHPropNewFunc     = MHPropNewFunc,
                     MHBlocks          = MHBlocks,
                     MHBlockNTimes     = MHBlockNTimes,
                     moveProbsList     = list(MH = 1.0),
                     moveNTimesList    = list(MH = 1),
                     nThin             = nThin,
                     saveFitness       = saveFitness,
                     verboseLevel      = verboseLevel,
                     ...)

    res$startingVal <- startingVal
    procFinal1(res)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MetropolisHastings <-
    function (nIters,              
              startingVal,
              logTarDensFunc,
              propNewFunc,
              logPropDensFunc,
              MHBlocks      = NULL,
              MHBlockNTimes = NULL,
              nThin         = 1,
              saveFitness   = FALSE, 
              verboseLevel  = 0,
              ...)    
{
    if (is.null(propNewFunc)) {
        MHPropNewFunc <- NULL
    }
    else {
        propNewFunc <- .check.func.do(propNewFunc, 
                                      argsReq   = c('block', 'currentDraw', '...'),
                                      retObjMsg = 'a numeric vector')
        MHPropNewFunc <-
            function (temperature, block, currentDraw, ...)
                propNewFunc(block, currentDraw, ...)
    }    

    if (is.null(logPropDensFunc)) {
        logMHPropDensFunc <- NULL
    }
    else {
        argsReq   <- c('block', 'currentDraw', 'proposalDraw', '...')
        retObjMsg <- 'a numeric value'
        logPropDensFunc <- .check.func.do(logPropDensFunc,
                                          argsReq   = argsReq,
                                          retObjMsg = retObjMsg)
        logMHPropDensFunc <-
            function (temperature, block, currentDraw, proposalDraw, ...)
                logPropDensFunc(block, currentDraw, proposalDraw, ...)
    }

    res <- TOEMCMain(nIters            = nIters,              
                     temperLadder      = c(1),
                     startingVals      = startingVal,
                     logTarDensFunc    = logTarDensFunc,
                     MHPropNewFunc     = MHPropNewFunc,
                     logMHPropDensFunc = logMHPropDensFunc,
                     MHBlocks          = MHBlocks,
                     MHBlockNTimes     = MHBlockNTimes,
                     moveProbsList     = list(MH = 1.0),
                     moveNTimesList    = list(MH = 1),
                     nThin             = nThin,
                     saveFitness       = saveFitness,
                     verboseLevel      = verboseLevel,
                     ...)

    res$startingVal <- startingVal
    procFinal1(res)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parallelTempering <-
    function (nIters,              
              temperLadder,
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc = NULL, 
              MHBlocks          = NULL,
              MHBlockNTimes     = NULL,
              moveProbsList     = NULL,
              moveNTimesList    = NULL,
              levelsSaveSampFor = NULL,
              nThin             = 1,
              saveFitness       = FALSE, 
              verboseLevel      = 0,
              ...)          
{
    if (is.null(moveProbsList))
        moveProbsList <- list(MH = 1.0)
    
    if (is.null(moveNTimesList)) {
        temperLadder   <- .check.temperLadder(temperLadder)
        nLevels        <- length(temperLadder)
        moveNTimesList <- list(MH = 1,
                               RE = nLevels)
    }

    movenamesAllowed <- c('MH', 'RE')
    .check.moveList.do(moveProbsList, implementedMovenames  = movenamesAllowed)
    .check.moveList.do(moveNTimesList, implementedMovenames = movenamesAllowed)
                       
    res <- TOEMCMain(nIters            = nIters,              
                     temperLadder      = temperLadder,
                     startingVals      = startingVals,
                     logTarDensFunc    = logTarDensFunc,
                     MHPropNewFunc     = MHPropNewFunc,
                     MHBlocks          = MHBlocks,
                     MHBlockNTimes     = MHBlockNTimes,
                     logMHPropDensFunc = logMHPropDensFunc, 
                     moveProbsList     = moveProbsList,
                     moveNTimesList    = moveNTimesList,
                     levelsSaveSampFor = levelsSaveSampFor,
                     nThin             = nThin,
                     saveFitness       = saveFitness,
                     verboseLevel      = verboseLevel,
                     ...)
    res
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evolMonteCarlo <-
    function (nIters,              
              temperLadder,
              startingVals,
              logTarDensFunc,
              MHPropNewFunc,
              logMHPropDensFunc = NULL, 
              MHBlocks          = NULL,
              MHBlockNTimes     = NULL,  
              moveProbsList     = NULL,
              moveNTimesList    = NULL,
              SCRWMNTimes       = NULL, 
              SCRWMPropSD       = NULL, 
              levelsSaveSampFor = NULL,
              nThin             = 1,
              saveFitness       = FALSE,
              verboseLevel      = 0,
              ...)    
{    
    if (is.null(moveNTimesList)) {
        temperLadder   <- .check.temperLadder(temperLadder)
        nLevels        <- length(temperLadder)
        moveNTimesList <- list(MH = 1,
                               RE = nLevels)
    }

    res <- TOEMCMain(nIters            = nIters,              
                     temperLadder      = temperLadder,
                     startingVals      = startingVals,
                     logTarDensFunc    = logTarDensFunc,
                     MHPropNewFunc     = MHPropNewFunc,
                     MHBlocks          = MHBlocks,
                     MHBlockNTimes     = MHBlockNTimes,
                     logMHPropDensFunc = logMHPropDensFunc, 
                     moveProbsList     = moveProbsList,
                     moveNTimesList    = moveNTimesList,
                     SCRWMNTimes       = SCRWMNTimes,
                     SCRWMPropSD       = SCRWMPropSD,
                     levelsSaveSampFor = levelsSaveSampFor,
                     nThin             = nThin,
                     saveFitness       = saveFitness,
                     verboseLevel      = verboseLevel,
                     ...)
    res
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.EMC <-
    function (x, ...)
{
    if (!is.null(x$temperLadder)) {
        cat('\nThe temperature ladder:\n')
        print(x$temperLadder, ...)
    }
    cat('\nThe overall acceptance rate summary:\n')
    print(x$acceptRatios, ...)
    cat('\n')
}
