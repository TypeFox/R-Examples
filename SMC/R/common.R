
###  $Id: common.R,v 1.4 2008/02/04 19:57:55 goswami Exp $
###  
###  File:    common.R
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

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.get.dimSummPerPeriod <-
    function (nStreams,
              dimPerPeriod,
              summaryFunc,
              ...)
{
    if (is.null(summaryFunc))
        return(dimPerPeriod)

    currentStreams    <- matrix(0, nrow = nStreams, ncol = dimPerPeriod)
    currentLogWeights <- numeric(nStreams)

    tmp <- summaryFunc(currentPeriod     = 1,
                       currentStreams    = currentStreams,
                       currentLogWeights = currentLogWeights,
                       ...)
    if (!is.vector(tmp)) {
        msg <- paste('The provided summaryFunc does not return a vector. ',
                     'The returned object:', sep = '')
        fatal(msg, tmp)
    }
    
    as.integer(length(tmp))
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.nStreamsPreResamp <-
    function (nStreamsPreResamp, nStreams)
{
    if (is.null(nStreamsPreResamp))
        return(as.integer(nStreams))

    nStreamsPreResamp <- .check.numericWithLLim(nStreamsPreResamp, 0)
    if (nStreamsPreResamp < nStreams) {
        msg <- paste('Please provide a valid value for "nStreamsPreResamp" :: ',
                     'it should be >= nStreams: ', nStreams, '. ',
                     'The given value:', sep = '')
        fatal(msg, nStreamsPreResamp)
    }

    as.integer(nStreamsPreResamp)
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

.check.propagateFunc <-
    function (propagateFunc)
{
    argsReq   <- c('currentPeriod', 'nStreamsToGenerate', 'lag1Streams',
                   'lag1LogWeights', 'startingStreams', '...')
    retObjMsg <- paste('a list with the following components:\n',
                       '   currentStreams: a numeric matrix\n',
                       'currentLogWeights: a numeric vector', sep = '')    
    .check.func.do(propagateFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.resampCriterionFunc <-
    function (resampCriterionFunc)
{
    if (is.null(resampCriterionFunc))
        return(NULL)

    argsReq   <- c('currentPeriod', 'currentStreams', 'currentLogWeights', '...')
    retObjMsg <- 'a logical value'    
    .check.func.do(resampCriterionFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.resampFunc <-
    function (resampFunc)
{
    if (is.null(resampFunc))
        return(NULL)

    argsReq   <- c('currentPeriod', 'currentStreams', 'currentLogWeights', '...')
    retObjMsg <- paste('a list with the following components:\n',
                       '   currentStreams: a numeric matrix\n',
                       'currentLogWeights: a numeric vector', sep = '')    
    .check.func.do(resampFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.summaryFunc <-
    function (summaryFunc)
{
    if (is.null(summaryFunc))
        return(NULL)

    argsReq   <- c('currentPeriod', 'currentStreams', 'currentLogWeights', '...')
    retObjMsg <- 'a numeric vector'
    .check.func.do(summaryFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.MHUpdateFunc <-
    function (MHUpdateFunc)
{
    argsReq   <- c('currentPeriod', 'nMHSteps', 'currentStreams', 'lag1Streams', 
                   'lag1LogWeights', '...')
    retObjMsg <- paste('a list with the following components:\n',
                       ' currentStreams: a numeric matrix\n',
                       'acceptanceRates: a numeric vector', sep = '')    
    .check.func.do(MHUpdateFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.generateStreamRepsFunc <-
    function (generateStreamRepsFunc)
{
    argsReq   <- c('currentPeriod', 'lag1Streams', 'lag1LogWeights',
                   'streamIndices', '...')
    retObjMsg <- 'a numeric matrix'
    .check.func.do(generateStreamRepsFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.generateNextStreamsFunc.PF <-
    function (generateNextStreamsFunc)
{
    argsReq   <- c('currentPeriod', 'lag1Streams', 'lag1LogWeights',
                   'streamIndices', 'startingStreams', '...')
    retObjMsg <- 'a numeric matrix'
    .check.func.do(generateNextStreamsFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.generateNextStreamsFunc.APF <-
    function (generateNextStreamsFunc)
{
    argsReq   <- c('currentPeriod', 'lag1Streams', 'lag1LogWeights',
                   'streamIndices', 'streamReps', 'startingStreams', '...')
    retObjMsg <- 'a numeric matrix'
    .check.func.do(generateNextStreamsFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

.check.logObsDensFunc <-
    function (logObsDensFunc)
{
    argsReq   <- c('currentPeriod', 'currentStreams', '...')
    retObjMsg <- 'a numeric vector'
    .check.func.do(logObsDensFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.SMCObj <-
    function (SMCObj)
{
    if (!inherits(SMCObj, 'SMC'))
        fatal('Please provide an object of class "SMC" for SMCObj')

    if (is.null(SMCObj$draws$streams)) {
        msg <- paste('Please re-run your SMC sampler with "returnStreams = TRUE"',
                     'and then pass the resulting SMC output object')
        fatal(msg)
    }
    
    if (is.null(SMCObj$draws$logWeights)) {
        msg <- paste('Please re-run your SMC sampler with "returnLogWeights = TRUE"',
                     'and then pass the resulting SMC output object')
        fatal(msg)
    }    
}

.check.generateObsFunc <-
    function (generateObsFunc)
{
    argsReq   <- c('currentPeriod', 'currentStreams', 'nPeriodsPast',
                   'lag1Obs', '...')
    retObjMsg <- 'a numeric vector'
    .check.func.do(generateObsFunc, argsReq = argsReq, retObjMsg = retObjMsg)    
}

.check.nPeriodsPast <-
    function (nPeriodsPast, nPeriods)
{
    if (nPeriodsPast > nPeriods) {
        msg <- paste('Please provide a valid nPeriodsPast :: it should be an ',
                     'integer <= nPeriods (= ', nPeriods, ')', sep = '')
        fatal(msg, nPeriodsPast)
    }

    as.integer(nPeriodsPast)
}

.check.weightedSummaryFunc <-
    function (weightedSummaryFunc)
{
    argsReq <- c('vals', 'weights', '...')
    retObjMsg <- 'a numeric value'
    .check.func.do(weightedSummaryFunc, argsReq = argsReq, retObjMsg = retObjMsg)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doCall <-
    function (func, argsList, dotsList)
{
    do.call(what = func, args = c(argsList, dotsList))
}
