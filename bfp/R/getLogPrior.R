#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getLogPrior.R] by DSB Mon 17/01/2011 15:37 (CET)>
##
## Description:
## Extract log prior from a model.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 29/11/2008   update for new model prior option
## 14/01/2011   update for model prior "dependent"
#####################################################################################

getLogPrior <- function (x       # a valid BayesMfp-Object of length 1 (otherwise only first element recognized)
                         )
{
    ## only take the first model
    x <- x[1]
    
    ## extract components
    ats <- attributes(x)   
    mp <- ats$priorSpecs$modelPrior
    
    powers <- x[[1]]$powers
    ucTerms <- x[[1]]$ucTerms
    
    ## depending on prior type, calculate the log prior value.
    if(identical(mp, "sparse"))
    {        
        maxs <- ats$shiftScaleMax[, 3]
        cards <- ats$shiftScaleMax[, 4]

        logVals <- numeric(length=length(powers))
        
        for (i in seq_along (powers))
        {
            deg <- length (powers[[i]])
            logVals[i] <- - lchoose(cards[i] - 1 + deg, deg) - log1p(maxs[i])
        }

        return (sum (logVals) - max (ats$indices$uc) * log (2))
    }
    else if(identical(mp, "dependent")) ## new (as of January 2011)
    {
        ## determine number of all covariates (covariate groups):
        nCovs <- length(ats$indices$ucList) + length(ats$indices$bfp)

        ## determine number of included covariates (covariate groups):
        fpDegrees <- sapply(powers, length)
        isIncludedFp <- fpDegrees > 0L
        
        nInclContinuous <- sum(isIncludedFp)
        nInclDiscrete <- length(ucTerms)
        
        nIncluded <- nInclContinuous + nInclDiscrete

        ## determine number of nonlinear covariates:
        isLinearFp <- sapply(powers,
                             identical,
                             y=1)
        isNonlinearFp <- isIncludedFp & (! isLinearFp)
        
        nNonlinear <- sum(isNonlinearFp)

        ## and determine number of possible nonlinear transformations
        ## for each continuous covariate
        nPossibleNonlinearTransforms <-
            if(nNonlinear == 0L)
                numeric(0)
            else
                sapply(ats$shiftScaleMax[isNonlinearFp, "maxDegree"],
                       getNumberPossibleFps) - 2L # subtract degree 0
                                        # and linear degree 1
        
        ## so the final result is:
        return(- log1p(nCovs) - lchoose(nCovs, nIncluded) -
               log1p(nInclContinuous) - lchoose(nInclContinuous, nNonlinear) -
               sum(log(nPossibleNonlinearTransforms)))
    }
    else ## flat
    {
        return(- max (ats$indices$uc) * log (2))
    }
}
