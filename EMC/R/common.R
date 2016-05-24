
###  $Id: common.R,v 1.17 2008/07/06 03:09:12 goswami Exp $
###  
###  File:    common.R
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
### Some utility functions 

fatal <-
    function (msg, var, formatMsg = TRUE,
              stopMsg = 'Fix the above problem first')
{
    if (formatMsg)
        cat(strsplit(msg, ' ')[[1]], '\n', fill = TRUE)
    else
        cat(msg, '\n')
    if (!missing(var))
        print(var)
    stop(stopMsg, call. = FALSE)
}

collectVarnames <-
    function (varnames, env = sys.frame(-1), simplify = FALSE, USE.NAMES = TRUE)
{
    sapply(varnames, function (vv) get(vv, envir = env),
           simplify = simplify, USE.NAMES = USE.NAMES)
}

procAcceptRatios1 <-
    function (samplerObj)
{
    aR <- as.data.frame(samplerObj$acceptRatios)
    rownames(aR) <- 'overall'
    aR[c('ratio', 'accepted', 'proposed')]
}

procAcceptRatios2 <-
    function (samplerObj)
{
    dAR <- samplerObj$detailedAcceptRatios$MH
    row.names(dAR) <- 'acceptRatios'
    dAR
}

procFinal1 <-
    function (res)
{
    res$levelsSaveSampFor    <- NULL
    res$temperLadder         <- NULL
    res$startingVals         <- NULL
    res$moveProbsList        <- NULL
    res$moveNTimesList       <- NULL

    res$draws                <- drop(res$draws)
    res$acceptRatios         <- procAcceptRatios1(res)
    res$detailedAcceptRatios <- procAcceptRatios2(res)
    res
}
   
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Some helper functions for TOEMCMain

.check.temperLadder <-
    function (temperLadder)
{
    msg <- paste('Please provide a valid temperLadder :: it should be a ',
                 'decreasing sequence of positive numbers. ',
                 'The given object:', sep = '')
    if (!is.vector(temperLadder) ||
        (length(temperLadder) == 0) ||
        !all(temperLadder > 0))
        fatal(msg, temperLadder)

    ord <- order(temperLadder, decreasing = TRUE)
    if (!identical(ord, seq_along(temperLadder))) 
        fatal(msg, temperLadder)
    
    as.double(temperLadder)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.startingValsVec <-
    function (startingVals, nLevels)
{
    if (!is.vector(startingVals))
        return(startingVals)

    if (length(startingVals) == 0) {
        msg <- paste('Please provide a vector for startingVals of length ',
                     'same as the dimension of the sample space. ',
                     'The given object:', sep = '')
        fatal(msg, startingVals)
    }
    
    matrix(rep(startingVals, nLevels), nLevels, byrow = TRUE)
}

.check.startingValsMat <-
    function (startingVals, nLevels)
{
    msg <- paste('Please provide a valid matrix for startingVals :: ',
                 'it should be a matrix with nLevels = ', nLevels, ' rows, ',
                 'where nLevels = length of the temperLadder. ',
                 'The given object:', sep = '')

    if (!is.matrix(startingVals))
        fatal(msg, startingVals)

    if ((nRows <- nrow(startingVals)) != nLevels) 
        fatal(msg, startingVals)    
    
    apply(startingVals, c(1, 2), as.double)
}

.check.startingVals <-
    function (startingVals, nLevels)
{
    startingVals <- .check.startingValsVec(startingVals, nLevels)
    .check.startingValsMat(startingVals, nLevels)
}
    
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.func.do <-
    function (func,
              argsReq,
              retObjMsg,
              funcname   = substitute(func),
              checkNames = TRUE)
{
    msg <- paste('Please provide a valid R function for "', funcname,
                 '" with argument(s):\n',
                 '(', toString(argsReq), ')\n\n',
                 'The function should return:\n',
                 retObjMsg, '\n\n',
                 'See the relevant help page for details.\n',
                 sep = '')

    if (!is.function(func)) {
        msg <- paste(msg, 'The given object:\n', sep = '')
        fatal(msg, func, formatMsg = FALSE)
    }

    argsHave <- names(formals(func))
    if (checkNames) {
        if (length(argsMissing <- setdiff(argsReq, argsHave)) > 0) {
            msg <- paste(msg, 'The following arguments are missing:', sep = '')
            fatal(msg, argsMissing, formatMsg = FALSE)
        }
        
        if (!('...' %in% argsReq) &&
            (length(argsExtra <- setdiff(argsHave, argsReq)) > 0)) {
            msg <- paste(msg, 'The function has the following extra arguments:',
                         sep = '')
            fatal(msg, argsExtra, formatMsg = FALSE)
        }
    }
    else {
        if (!('...' %in% argsReq) &&
            ((nHave <- length(argsHave)) != (nReq <- length(argsReq)))) {
            msg <- paste(msg, 'The given function takes ', nHave, ' arguments ',
                         'instead of ', nReq, sep = '')
            fatal(msg, formatMsg = FALSE)
        }
    }
    
    as.function(func)
}

.check.logTarDensFunc <-
    function (logTarDensFunc)
{
    argsReq   <- c('draw', '...')
    retObjMsg <- 'a numeric value'        
    .check.func.do(logTarDensFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.MHPropNewFunc <-
    function (MHPropNewFunc)
{
    argsReq   <- c('temperature', 'block', 'currentDraw', '...')
    retObjMsg <- 'a numeric vector'
    .check.func.do(MHPropNewFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.logMHPropDensFunc <-
    function (logMHPropDensFunc)
{
    if (is.null(logMHPropDensFunc))
        return(NULL)

    argsReq   <- c('temperature', 'block', 'currentDraw', 'proposalDraw', '...')
    retObjMsg <- 'a numeric value'
    .check.func.do(logMHPropDensFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.MHRWMPropDispVec <-
    function (MHRWMPropDisp, sampDim)
{
    if (!is.vector(MHRWMPropDisp))
        return(MHRWMPropDisp)
    
    if (length(MHRWMPropDisp) != sampDim) {
        msg <- paste('Please provide a vector for MHRWMPropDisp :: ',
                     'it should be a vector of length:\n',
                     'sampDim  =  ', sampDim, ', where\n',
                     'sampDim = number of columns of the starting value matrix. ',
                     'The given vector:', sep = '')
        fatal(msg, MHRWMPropDisp, formatMsg = FALSE)
    }
    
    diag(MHRWMPropDisp)
}

.check.MHRWMPropDispMat <-
    function (MHRWMPropDisp, sampDim, nLevels, temperLadder)
{
    if (!is.matrix(MHRWMPropDisp))
        return(MHRWMPropDisp)

    if ((dimHave <- dim(MHRWMPropDisp)) != c(sampDim, sampDim)) {
        msg <- paste('Please provide a matrix for MHRWMPropDisp :: ',
                     'it should be a matrix of dimension:\n',
                     'sampDim x sampDim  =  ',
                     sampDim, ' x ', sampDim, ', where\n',
                     'sampDim = number of columns of the starting value matrix. ',
                     'The given matrix is of dimension:', sep = '')
        fatal(msg, dimHave, formatMsg = FALSE)
    }
    
    arr <- array(dim = c(sampDim, sampDim, nLevels))
    for (ii in seq_len(nLevels))
        arr[ , , ii] <- temperLadder[ii] * MHRWMPropDisp
    
    arr
}

.check.MHRWMPropDispArr <-
    function (MHRWMPropDisp, sampDim, nLevels)
{
    msg <- paste('Please provide an array of matrices for MHRWMPropDisp ::\n',
                 'it should be an array of matrices of dimension:\n',
                 'nLevels x sampDim x sampDim  =  ',
                 nLevels, ' x ', sampDim, ' x ', sampDim, ', where\n',
                 'nLevels = length of the temperLadder, and\n',
                 'sampDim = number of columns of the starting value matrix.\n',
                 sep = '')

    if (!is.array(MHRWMPropDisp)) 
        fatal(msg, MHRWMPropDisp, formatMsg = FALSE)

    if ((dimHave <- dim(MHRWMPropDisp)) != c(sampDim, sampDim, nLevels)) {    
        msg <- paste(msg, 'The given array is of dimension:', sep = '')
        fatal(msg, dimHave, formatMsg = FALSE)
    }

    apply(MHRWMPropDisp, c(1, 2, 3), as.double)
}

.check.MHRWMPropDisp <-
    function (MHRWMPropDisp, sampDim, nLevels, temperLadder, MHPropNewFunc)
{
    if (is.null(MHRWMPropDisp)) {
        if (is.null(MHPropNewFunc)) {
            msg <- paste('Please either provide a valid MHPropNewFunc or\n',
                         'an array of matrices for MHRWMPropDisp ::\n',
                         'it should be an array of matrices of dimension:\n',
                         'nLevels x sampDim x sampDim  =  ',
                         nLevels, ' x ', sampDim, ' x ', sampDim, ', where\n',
                         'nLevels = length of the temperLadder, and\n',
                         'sampDim = number of columns of the starting value matrix.\n',
                         sep = '')            
            fatal(msg, MHRWMPropDisp, formatMsg = FALSE)
        }
        return(NULL)
    }

    MHRWMPropDisp <- .check.MHRWMPropDispVec(MHRWMPropDisp, sampDim = sampDim)
    MHRWMPropDisp <- .check.MHRWMPropDispMat(MHRWMPropDisp,
                                             sampDim      = sampDim,
                                             nLevels      = nLevels,
                                             temperLadder = temperLadder)
    .check.MHRWMPropDispArr(MHRWMPropDisp,
                            sampDim = sampDim,
                            nLevels = nLevels)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.MHBlocks <-
    function (MHBlocks, sampDim)
{
    if (is.null(MHBlocks)) 
        return(lapply(seq_len(sampDim), as.integer))

    if (!is.list(MHBlocks) ||
        (length(MHBlocks) == 0)) {
        msg <- paste('Please provide a list for MHBlocks :: it should be a list',
                     'of blocks of dimensions to be sampled together.',
                     'The following was provided instead:')
        fatal(msg, MHBlocks)
    }
    
    dimsGiven <- unlist(MHBlocks)
    dimsNotGiven <- setdiff(seq_len(sampDim), dimsGiven)
    if (length(dimsNotGiven) > 0) {
        msg <- paste('Please provide a list for MHBlocks :: it should be a list',
                     'of blocks, the union of which should cover all the',
                     'dimensions of the sample space. The following dimensions',
                     'were not covered by the MHBlocks:')
        fatal(msg, MHBlocks)
    }
    
    lapply(MHBlocks, as.integer)
}

.check.MHBlockNTimes <-
    function (MHBlockNTimes, nMHBlocks)
{
    if (is.null(MHBlockNTimes))
        return(as.integer(rep(1, nMHBlocks)))

    if (!is.vector(MHBlockNTimes) ||
        (length(MHBlockNTimes) != nMHBlocks) ||
        (any(MHBlockNTimes <= 0))) {
        msg <- paste('Please provide a vector for MHBlockNTimes :: it should ',
                     'be of length nMHBlocks = ', nMHBlocks, ', where ',
                     'nMHBlocks = length of MHBlocks list. The following was ',
                     'provided instead:')
        fatal(msg, MHBlockNTimes)
    }
    
    as.integer(MHBlockNTimes)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.isMHRWMVec <-
    function (isMHRWMVec, nMHBlocks)
{
    if (is.null(isMHRWMVec))
        return(NULL)

    if (!all(is.logical(isMHRWMVec)) ||
        (length(isMHRWMVec) != nMHBlocks)) {
        msg <- paste('Please provide a vector for isMHRWMVec :: it should ',
                     'be of length nMHBlocks = ', nMHBlocks, ', where ',
                     'nMHBlocks = length of MHBlocks list. The following was ',
                     'provided instead:')
        fatal(msg, isMHRWMVec)
    }
    
    as.logical(isMHRWMVec)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.get.implementedMovenames <-
    function ( )
    c('MH', 'RC', 'SC',
      'RE', 'BCE', 'BIRE', 'BSE', 'CE') 

.get.specialMovenames <-
    function ( )
    c('BCE', 'BIRE', 'BSE', 'CE') 

.get.crossoverMovenames <-
    function ( )
    c('RC', 'SC')

.get.implementedMovesList <-
    function (func, implementedMovenames = .get.implementedMovenames( ))
    sapply(implementedMovenames, function (cc) func(0), simplify = FALSE)

.check.moveList.do <-
    function (aList,
              listname             = substitute(aList),
              implementedMovenames = .get.implementedMovenames( ))
{
    msg <- paste('Please provide a list for', listname)

    if (!is.list(aList)) {
        msg <- paste(msg, ', given object:', sep = '')
        fatal(msg, aList)
    }

    given <- names(aList)
    if (!('MH' %in% given)) {
        msg <- paste(msg, 'having a component called MH, given list:')
        fatal(msg, aList)
    }

    if (length(notKnown <- setdiff(given, implementedMovenames)) > 0) {
        msg <- paste('Some of the components of the list', listname, 'are not',
                     'implemented move names, allowable names are:')
        fatal(msg, implementedMovenames)
    }

    if (any(unlist(aList) < 0)) {
        msg <- paste('All the components of', listname, 'are not positive,',
                     'given list:')
        fatal(msg, aList)
    }

    if (aList$MH <= 0) {
        msg <- paste('The "MH" component of', listname, 'should be positive,',
                     'given list:')
        fatal(msg, aList)
    }
}

.check.moveList.fill <-
    function (aList, func)
{
    ll <- .get.implementedMovesList(func)
    for (nn in names(aList))
        ll[[nn]] <- (func)(aList[[nn]])

    ll
}    

.check.moveProbsList <-
    function (moveProbsList, nLevels)
{
    if (is.null(moveProbsList)) {
        if (nLevels == 1) {
            moveProbsList <- list(MH = as.double(1))
        }
        else {
            moveProbsList <- list(MH = as.double(0.4),
                                  RC = as.double(0.3),
                                  SC = as.double(0.3))
        }
        return(.check.moveList.fill(moveProbsList, as.double))
    }

    .check.moveList.do(moveProbsList)
    if (abs(sum(unlist(moveProbsList)) - 1) > (getOption('ts.eps') / 2)) {
        msg <- paste('The probabilities should add up to 1 in moveProbsList,',
                     'given list:')
        fatal(msg, moveProbsList)
    }

    .check.moveList.fill(moveProbsList, as.double)
}

.check.moveNTimesList <-
    function (moveNTimesList, moveProbsList, nLevels)
{
    nmpl     <- names(moveProbsList[moveProbsList > 0])
    defaults <- list(MH   = as.integer(1),
                     SC   = as.integer(floor(nLevels / 2)),
                     RC   = as.integer(floor(nLevels / 2)),
                     RE   = as.integer(nLevels),
                     BCE  = as.integer(max(1, nLevels - 3)),
                     BIRE = as.integer(max(1, nLevels - 3)),
                     BSE  = as.integer(max(1, nLevels - 3)),
                     CE   = as.integer(max(1, nLevels - 3)))
    
    if (is.null(moveNTimesList)) {
        if (nLevels == 1) {
            moveNTimesList <- list(MH = as.integer(1))
        }
        else {
            moveNTimesList <- defaults[nmpl] 
        }
        return(.check.moveList.fill(moveNTimesList, as.integer))
    }

    .check.moveList.do(moveNTimesList)
    ll <- .check.moveList.fill(moveNTimesList, as.integer)

    nmntl <- names(moveNTimesList[moveNTimesList > 0])    
    if (length(forgot <- setdiff(nmpl, nmntl)) > 0) {
        warning('Some of the components (',  toString(forgot), ') of ',
                'moveNTimesList were made positive since the corresponding ',
                'components of moveProbsList were positive', call. = FALSE)
        for (nn in forgot)
            ll[[nn]] <- defaults[[nn]]
    }
    
    for (nn in intersect(names(moveNTimesList), .get.specialMovenames( ))) {
        if (moveNTimesList[[nn]] > nLevels - 2) {
            msg <- paste('Please provide a valid moveNTimes for move ', nn,
                         ' in moveNTimesList :: it should have non-negative',
                         ' number <= ', nLevels - 2, ' since nLevels = ',
                         nLevels, '. The given moveNTimes:', sep = '')
            fatal(msg, moveNTimesList[[nn]])
        }
    }

    ll
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.SCControls.do <-
    function (varname, var, default, moveNTimesList)
{
    if (moveNTimesList$SC == 0)
        return(0)

    if (is.null(var))
        var <- default
    if (is.null(var) || (var < 0)) {
        msg <- paste('Please provide a positive value for ', varname,
                     '. The given object:', sep = '')
        fatal(msg, var)
    }
    
    var
}           

.check.SCRWMPropSD <-
    function (SCRWMPropSD, moveNTimesList)
{
    as.double(.check.SCControls.do(varname        = 'SCRWMPropSD',
                                   var            = SCRWMPropSD,
                                   default        = NULL,
                                   moveNTimesList = moveNTimesList))
}

.check.SCRWMNTimes <-
    function (SCRWMNTimes, moveNTimesList)
{
    as.integer(.check.SCControls.do(varname        = 'SCRWMNTimes',
                                    var            = SCRWMNTimes,
                                    default        = 1,
                                    moveNTimesList = moveNTimesList))
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.get.selectionCodes <-
    function ( )
    c('random', 'best', 'worst')

.get.defaults.selectionCodes <-
    function ( )
{
    list(RC = c('best', 'best'),
         SC = c('best'))
}

.check.mSCList.do <-
    function (moveSelectionCodesList, moveNTimesList, movename, nVals,
              selectionCodes = .get.selectionCodes( ))
{
    vv   <- moveSelectionCodesList[[movename]]
    nn   <- moveNTimesList[[movename]]
    msg1 <- paste('The', movename, 'component of moveSelectionCodesList ')
    
    if (is.null(vv)) {
        if (nn > 0) {
            msg <- paste(msg1, 'is not provided but the corresponding component ',
                         'of moveNTimesList is positive: ', nn)
            fatal(msg)
        }

        return(invisible(0))
    }
    
    if (length(vv) != nVals) {
        msg <- paste(msg1, ' should be of length ', nVals, ', given vector:',
                     sep = '')
        fatal(msg, vv)
    }

    if (length(invalid <- setdiff(tolower(vv), selectionCodes)) > 0) {
        msg <- paste(msg1, ' has some invalid selection codes: ',
                     toString(invalid), ', allowed codes are:', sep = '')
        fatal(msg, selectionCodes)
    }
}

.check.mSCList <-
    function (moveSelectionCodesList, moveNTimesList)
{
    defaults <- .get.defaults.selectionCodes( )

    if (is.null(moveSelectionCodesList))
        return(sapply(defaults, as.character, simplify = FALSE))

    if (!is.list(moveSelectionCodesList)) {
        msg <- paste('Please provide a list for moveSelectionCodesList,',
                     'given object:')
        fatal(msg, moveSelectionCodesList)
    }
    
    for (vv in names(defaults)) {
        .check.mSCList.do(moveSelectionCodesList = moveSelectionCodesList,
                          moveNTimesList         = moveNTimesList,
                          movename               = vv,
                          nVals                  = length(defaults[[vv]]))
    }
    
    needComps  <- names(defaults)
    nmscl      <- names(moveSelectionCodesList)
    
    if (length(extraComps <- setdiff(nmscl, needComps)) > 0) 
        warning('The following extra components of moveSelectionCodesList ',
                'will be ignored: ', toString(extraComps), call. = FALSE)

    missingComps           <- setdiff(needComps, nmscl)
    moveSelectionCodesList <- c(moveSelectionCodesList, defaults[missingComps])
    sapply(moveSelectionCodesList[needComps], as.character, simplify = FALSE)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.get.defaults.selectionTempers <-
    function (temperColdest)
{
    list(RC  = rep(temperColdest, 2),
         SC  = temperColdest,
         BCE = temperColdest)
}

.check.mSTList.do <-
    function (moveSelectionTempersList, moveNTimesList, movename, nVals)
{
    vv   <- moveSelectionTempersList[[movename]]
    nn   <- moveNTimesList[[movename]]
    msg1 <- paste('The', movename, 'component of moveSelectionTempersList')

    if (is.null(vv)) {
        if (nn > 0) {
            msg <- paste(msg1, 'is not provided but the corresponding component ',
                         'of moveNTimesList is positive: ', nn)
            fatal(msg)
        }
        
        return(invisible(0))
    }
    
    if (length(vv) != nVals) {
        msg <- paste(msg1, ' should be of length ', nVals, ', given vector:',
                     sep = '')
        fatal(msg, vv)
    }
    
    if (any(bad <- (vv < 0))) {
        msg <- paste(msg1, 'has some non-negative selectionTempers:')
        fatal(msg, vv[bad])
    }
}

.check.mSTList <-
    function (moveSelectionTempersList, moveNTimesList, temperColdest)
{
    defaults <- .get.defaults.selectionTempers(temperColdest)
    
    if (is.null(moveSelectionTempersList))
        return(sapply(defaults, as.double, simplify = FALSE))

    if (!is.list(moveSelectionTempersList)) {
        msg <- paste('Please provide a list for moveSelectionTempersList,',
                     'given object:')
        fatal(msg, moveSelectionTempersList)
    }

    for (vv in names(defaults)) {
        .check.mSTList.do(moveSelectionTempersList = moveSelectionTempersList,
                          moveNTimesList           = moveNTimesList,
                          movename                 = vv,
                          nVals                    = length(defaults[[vv]]))
    }    

    needComps <- names(defaults)
    nmstl     <- names(moveSelectionTempersList)
    
    if (length(extraComps <- setdiff(nmstl, needComps)) > 0)
        warning('The following extra components of moveSelectionTempersList ',
                'will be ignored: ', toString(extraComps), call. = FALSE)

    missingComps             <- setdiff(needComps, nmstl)
    moveSelectionTempersList <- c(moveSelectionTempersList, defaults[missingComps])
    sapply(moveSelectionTempersList[needComps], as.double, simplify = FALSE)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.levelsSaveSamplesFor <-
    function (levelsSaveSampFor, nLevels)
{
    if (is.null(levelsSaveSampFor))
        return(as.integer(nLevels))

    if (!is.vector(levelsSaveSampFor) ||
        (length(levelsSaveSampFor) > nLevels) ||
        !all((1 <= levelsSaveSampFor) &&
             (levelsSaveSampFor <= nLevels))) {
        msg <- paste('Please provide an integer vector for levelsSaveSampFor :: ',
                     'it should be of length  <=  nLevels (= ', nLevels, '), ',
                     'where nLevels = length of the temperLadder.', sep = '')
        fatal(msg, levelsSaveSampFor)
    }
    
    as.integer(levelsSaveSampFor)    
}

.check.logical <-
    function (val, varname = substitute(val))
{
    if (!is.logical(val)) {
        msg <- paste('Please provide a logical value for "', varname, '". ',
                     'The given object:', sep = '') 
        fatal(msg, val)
    }
    
    as.logical(val)
}

.check.numericWithLLim <-
    function (val, llim, varname = substitute(val), retFunc = as.integer)
{
    if (!is.numeric(val)) {        
        msg <- paste('Please provide an integer value for "', varname, '". ',
                     'The given object:', sep = '')
        fatal(msg, val)
    }

    if (!is.na(llim) &&
        (val < llim)) {
        msg <- paste('Provided "', varname, '" is too low :: it should be ',
                     'an integer >= ', llim, '. The given value:', sep = '')
        fatal(msg, val)
    }

    retFunc(val)
}

.check.timeInSecs <-
    function (timeInSecs)
{
    timeInSecs <- .check.numericWithLLim(timeInSecs, NA, retFunc = as.numeric)

    if (timeInSecs > 0)
        warning('Since timeInSecs: ', timeInSecs, ', which is positive, ',
                'nIters will be ignored', call. = FALSE)
    
    as.double(timeInSecs)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.statsFuncList <-
    function (statsFuncList)
{
    if (!is.list(statsFuncList)) {
        msg <- paste('Please provide a list of functions of one argument',
                     '(xx) each for statsFuncList, given object:')
        fatal(msg, statsFuncList)
    }
    
    argsReq   <- 'xx'
    retObjMsg <- 'a numeric value'
    for (ii in seq_along(statsFuncList))
        statsFuncList[[ii]] <-
            .check.func.do(statsFuncList[[ii]],
                           argsReq   = argsReq,
                           retObjMsg = retObjMsg,
                           funcname  = paste('statsFuncList[[', ii, ']]', sep = ''))

    statsFuncList
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buildLadder <-
    function (limits,
              length,
              scheme,
              param,
              plotit = FALSE)
{
    sl   <- sort(limits, decreasing = TRUE)
    ulim <- sl[1]
    llim <- sl[2]

    implementedSchemes <- c('linear',
                            'log',
                            'geometric',
                            'mult-power',
                            'add-power',                              
                            'reciprocal',
                            'exponential',
                            'tangent')
    scheme             <- tolower(scheme)
    scheme             <- agrep(scheme, implementedSchemes, value = TRUE)
    if (length(scheme) != 1)
        stop('Please provied a valid ladder builder scheme, ',
             'implemented schemes are:\n',
             paste('linear ',
                   'log ',
                   'geometric ',
                   'mult-power ',
                   'add-power   (needs a positive param) ',
                   'reciprocal ',
                   'exponential (needs a positive param) ',
                   'tangent     (needs a positive param)',
                   sep = '\n'))

    specialSchemes <- c('exponential', 'tangent')
    if ((llim < 1.0) &&
        (scheme %in% specialSchemes))
        stop('Please choose a different scheme for building the temperature ',
             'ladder ::\n the lower temperature limit cannot be < 1 for schemes ',
             'in [ ', toString(specialSchemes), ' ]')
    
    ladderLinear <- seq(ulim, llim, len = length)
    if (identical(scheme, 'linear')) {
        ladder <- ladderLinear
    }
    else if (identical(scheme, 'log')) {
        ladder <- exp(seq(log(ulim), log(llim), len = length))
    }
    else if (identical(scheme, 'geometric')) {
        ladder <- llim * exp((log(ulim / llim) / (length - 1)) *
                             seq.int(length - 1, 0))
    }
    else if (identical(scheme, 'mult-power')) {
        alpha  <- log(ulim / llim) / log(length)
        ladder <- llim * exp(alpha * seq.int(length, 1))
    }
    else if (identical(scheme, 'add-power')) {
        stopifnot(length(param) > 0)
        alpha  <- param[1]
        ladder <- llim + (((ulim - llim) / (length - 1)^alpha) *
                          (seq.int(length - 1, 0)^alpha))
    }
    else if (identical(scheme, 'reciprocal')) {
        ladder <- 1.0 / seq(1.0 / ulim, 1.0 / llim, len = length)
    }
    else if (identical(scheme, 'exponential')) {
        stopifnot(length(param) > 0)
        alpha  <- param[1]
        stopifnot(alpha > 0)
        ladder <- seq((log(ulim))^(1 / alpha),
                      (log(llim))^(1 / alpha), len = length)
        ladder <- exp(ladder^alpha)
    }
    else if (identical(scheme, 'tangent')) {
        stopifnot(length(param) > 0)
        alpha  <- param[1]
        ladder <- seq(atan(ulim^(1 / alpha)),
                      atan(llim^(1 / alpha)), len = length)
        ladder <- tan(ladder)^alpha
    }
    
    if (plotit) {
        plot(ladderLinear, ladder, xlab = 'linear', ylab = scheme)
        abline(a = 0.0, b = 1.0)
    }

    ladder
}

.get.temperLadder <-
    function (temperLadder,
              temperLimits,
              ladderLen,
              scheme,
              schemeParam)
{
    if (!is.null(temperLadder))
        return(.check.temperLadder(temperLadder))

    if (is.null(temperLimits) ||
        !is.vector(temperLimits) ||
        (length(temperLimits) != 2) ||
        !all(temperLimits > 0) ||
        (temperLimits[1] == temperLimits[2])) {
        msg <- paste('Please provide either a non-null temperLadder or',
                     'temperLimits with two distinct positive entries,',
                     'given limits:')
        fatal(msg, temperLimits)
    }

    ladderLen <- .check.numericWithLLim(ladderLen, 1)

    as.double(buildLadder(limits = temperLimits,
                          length = ladderLen,
                          scheme = scheme,
                          param  = schemeParam,
                          plotit = FALSE))    
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.acceptRatioLimits <-
    function (acceptRatioLimits)
{
    if (is.null(acceptRatioLimits) ||
        !is.vector(acceptRatioLimits) ||
        length(acceptRatioLimits) != 2 ||
        any(acceptRatioLimits < 0.0) ||
        any(acceptRatioLimits > 1.0)) {
        msg <- paste('Please provide a (increasing) vector for two elements for',
                     'acceptRatioLimits; given object:')
        fatal(acceptRatioLimits)
    }
    
    as.double(acceptRatioLimits)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Some functions for processing the acceptRatiosList

.listToDF <-
    function (ll)
{
    do.call(what = rbind, args = ll)
}

.get.detailedAcceptRatios.MH <-
    function (comp, roundBy)
{
    arl <- lapply(seq_len(nrow(comp)),
                  function (ii)
              {
                  rr <- comp[ii, ]
                  oo <- order(rr, decreasing = FALSE)
                  nn <- length(oo)
                  data.frame(argminBlock = oo[1], min = rr[oo[1]],
                             argmaxBlock = oo[nn], max = rr[oo[nn]])
              })
    
    arm           <- .listToDF(arl)
    rownames(arm) <- rownames(comp)
    list(MH = roundDF(arm, roundBy = roundBy))
}

.get.row.col <-
    function (index, nRows)
{
    rr <- index %% nRows
    cc <- index %/% nRows
    if (rr == 0) {
        rr <- nRows
    }
    else {
        cc <- cc + 1
    }
    c(rr, cc)
}

.get.min.max.mat <-
    function (mat)
{
    nRows <- nrow(mat)
    xx <- c(mat)
    oo <- order(xx, decreasing = FALSE, na.last = NA)
    nn <- length(oo)
    if (nn == 0) {
        return(list(argmin = NA,
                    min    = NA,
                    argmax = NA,
                    max    = NA))
    }
    
    argmin <- .get.row.col(oo[1], nRows)
    argmax <- .get.row.col(oo[nn], nRows)
    list(argmin = argmin,
         min    = xx[oo[1]],
         argmax = argmax,
         max    = xx[oo[nn]])    
}

.str1 <-
    function (ll)
{
    with(ll,
         data.frame(argminLevels = sprintf("%2d <--> %2d", argmin[1], argmin[2]),
                    min          = min,
                    argmaxLevels = sprintf("%2d <--> %2d", argmax[1], argmax[2]),
                    max          = max))
}

.get.detailedAcceptRatios.various <-
    function (comp)
{
    bal <- .get.min.max.mat(comp)

    targetLevel    <- nrow(comp)
    twol           <- .get.min.max.mat(comp[targetLevel, , drop = FALSE])
    twol$argmin[1] <- targetLevel
    twol$argmax[1] <- targetLevel

    list(BetweenAllLevels      = .str1(bal),
         TargetWithOtherLevels = .str1(twol))
}

roundDF <-
    function (adf, roundBy)
{
    for (vv in names(adf)) {
        if (!is.numeric(adf[[vv]]))
            next

        adf[[vv]] <- round(adf[[vv]], digits = roundBy)
    }
    adf
}

.get.detailedAcceptRatios.CE <-
    function (comp, roundBy)
{
    if (is.null(comp))
        return(NULL)

    arl <- lapply(seq_len(nrow(comp)),
                  function (ii)
              {
                  rr <- comp[ii, ]
                  oo <- order(rr, decreasing = FALSE, na.last = NA)
                  nn <- length(oo)
                  data.frame(argminSelLength = oo[1], min = rr[oo[1]],
                             argmaxSelLength = oo[nn], max = rr[oo[nn]])
              })
    
    arm           <- .listToDF(arl)
    rownames(arm) <- rownames(comp)
    list(CE = roundDF(arm, roundBy = roundBy))
}

.get.detailedAcceptRatios <-
    function (acceptRatiosList, roundBy = 3)
{
    variousList <- NULL
    for (vv in setdiff(names(acceptRatiosList), c('MH', 'CE'))) {
        tmp <- .get.detailedAcceptRatios.various(acceptRatiosList[[vv]])

        rownames(tmp$BetweenAllLevels)      <- vv
        rownames(tmp$TargetWithOtherLevels) <- vv
        
        if (is.null(variousList)) {
            variousList <- tmp
            next
        }
        
        variousList$BetweenAllLevels      <- rbind(variousList$BetweenAllLevels,
                                                   tmp$BetweenAllLevels)
        variousList$TargetWithOtherLevels <- rbind(variousList$TargetWithOtherLevels,
                                                   tmp$TargetWithOtherLevels)        
    }

    variousList$BetweenAllLevels      <- roundDF(variousList$BetweenAllLevels,
                                                 roundBy = roundBy)
    variousList$TargetWithOtherLevels <- roundDF(variousList$TargetWithOtherLevels,
                                                 roundBy = roundBy)

    c(.get.detailedAcceptRatios.MH(acceptRatiosList$MH, roundBy = roundBy),
      variousList,
      .get.detailedAcceptRatios.CE(acceptRatiosList$CE, roundBy = roundBy))
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doCall <-
    function (func, argsList, dotsList)
{
    do.call(what = func, args = c(argsList, dotsList))
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.oneIterMHFuncGen.checkLogDens <-
    function (temperature, block, draw, logDens, level, iter)
{
    tmp <- collectVarnames(c('iter', 'level', 'temperature', 'block'),
                           simplify = TRUE)
    cat('\n')
    print(data.frame(value = tmp))
    cat('The current draw:\n')
    print(draw)
    cat('The log density:\n')
    print(logDens)
    fatal(paste('In MH, alpha is non finite, or is not of length 1,',
                'inspect the above numbers'))
}


.oneIterMHFuncGen.checkPropNew <-
    function (temperature, block, curr, prop, level, iter)
{
    tmp <- collectVarnames(c('iter', 'level', 'temperature', 'block'),
                           simplify = TRUE)
    print(data.frame(value = tmp))
    cat('The current draw:\n')
    print(curr)
    cat('The proposal draw:\n')
    print(prop)
    fatal(paste('In MH, current and proposal draw are not of the same length,',
                'inspect the above numbers and the MHPropNewFunc'))
}


.oneIterMHFuncGen.checkPropDens <-
    function (temperature, block, curr, prop, logPropDens, level, iter)
{
    tmp <- collectVarnames(c('iter', 'level', 'temperature', 'block', 'logPropDens'),
                           simplify = TRUE)
    print(data.frame(value = tmp))
    cat('The current draw:\n')
    print(curr)
    cat('The proposal draw:\n')
    print(prop)
    cat('The log proposal density:\n')
    print(logPropDens)
    fatal(paste('In MH, logDens of proposal draw is not of length 1,',
                'inspect the above numbers and the logMHPropDensFunc'))
}


.oneIterMHFuncGen.checkAlpha <-
    function(temperature, block, curr, prop, level, iter,
             logDenominators, logNumerators, alpha)
{
    tmp <- collectVarnames(c('iter', 'level', 'temperature', 'block', 'logDens'),
                           simplify = TRUE)
    print(data.frame(value = tmp))
    cat('The current draw: logDens:', logDenominators[1], '\n')
    print(curr)
    cat('The proposal draw: logDens:', logNumerators[1], '\n')
    print(prop)
    fatal('In MH, alpha is non finite, inspect the above numbers')                
}


.oneIterMHFuncGen <-
    function (lTDF, MHPNF, lMHPDF)
{
    func <- 
        function (temperature, block, draw, logDens, level, iter, ...)
        {        
            logDenominators    <- logNumerators <- c(0, 0)
            logDenominators[1] <- logDens / temperature
            
            prop <- MHPNF(temperature = temperature,
                          block       = block,
                          currentDraw = draw, ...)
            if (length(draw) != length(prop))
                .oneIterMHFuncGen.checkPropNew(temperature = temperature,
                                               block       = block,
                                               curr        = draw,
                                               prop        = prop,
                                               level       = level,
                                               iter        = iter)

            logPropDens <- lTDF(prop, ...)
            if ((length(logPropDens) != 1))
                .oneIterMHFuncGen.checkLogDens(temperature = temperature,
                                               block       = block,
                                               draw        = draw,
                                               logDens     = logPropDens,
                                               level       = level,
                                               iter        = iter)                
            logNumerators[1] <- logPropDens / temperature
            
            if (!is.null(lMHPDF)) {
                logPropDens <- lMHPDF(temperature  = temperature,
                                      block        = block,
                                      currentDraw  = draw,
                                      proposalDraw = prop, ...)
                if ((length(logPropDens) != 1) || !is.finite(logPropDens))
                    .oneIterMHFuncGen.checkPropDens(temperature = temperature,
                                                    block       = block,
                                                    curr        = draw,
                                                    prop        = prop,
                                                    logPropDens = logPropDens,
                                                    level       = level,
                                                    iter        = iter)
                logDenominators[2] <- logPropDens
                
                logPropDens <- lMHPDF(temperature  = temperature,
                                      block        = block,
                                      currentDraw  = prop,
                                      proposalDraw = draw, ...)
                if ((length(logPropDens) != 1) || !is.finite(logPropDens))
                    .oneIterMHFuncGen.checkPropDens(temperature = temperature,
                                                    block       = block,
                                                    curr        = prop,
                                                    prop        = draw,
                                                    logPropDens = logPropDens,
                                                    level       = level,
                                                    iter        = iter)
                logNumerators[2] <- logPropDens
            }
            
            alpha <- min(1.0, exp(sum(logNumerators - logDenominators)))        
            if (!is.finite(alpha))
                .oneIterMHFuncGen.checkAlpha(temperature     = temperature,
                                             block           = block,
                                             curr            = draw,
                                             prop            = prop,
                                             level           = level,
                                             iter            = iter,
                                             logDenominators = logDenominators,
                                             logNumerators   = logNumerators,
                                             alpha           = alpha)
            list(prop = prop, logPropDens = logPropDens, alpha = alpha)
        }
    as.function(func)
}
