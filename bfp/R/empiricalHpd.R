#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[empiricalHpd.R] by DSB Mon 05/10/2009 10:22 (CEST)>
##
## Description:
## Compute a MC HPD estimate for a scalar parameter.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 01/10/2009   add argument checks, names for return vector
## 05/10/2009   comments, some beautifications
#####################################################################################

empiricalHpd <- function (theta,        # sample vector of parameter
                          level         # credible level
                          )
{
    ## check that theta is numeric, and that level is in (0, 1)
    stopifnot(is.numeric(theta),
              0 < level && 1 > level)

    ## how many samples are saved in theta?
    nSamples <- length (theta)

    ## get the sorted samples vector
    thetaSorted <- sort.int (theta, method = "quick")

    ## how many different credible intervals with "level"
    ## do we need to compare?
    nIntervals <- ceiling (nSamples * (1 - level))

    ## these are the start indexes of the intervals
    startIndexes <- seq_len (nIntervals)

    ## and these are the end indexes of the intervals
    endIndexes <- nSamples - nIntervals + startIndexes
    
    ## which interval has the smallest range?
    smallestInterval <- which.min(thetaSorted[endIndexes] - thetaSorted[startIndexes])

    ## then return the bounds of this smallest interval
    return (c (lower=thetaSorted[startIndexes[smallestInterval]],
               upper=thetaSorted[endIndexes[smallestInterval]]))
}
