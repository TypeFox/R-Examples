#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[scrHpd.R] by DSB Mit 26/01/2011 14:03 (CET)>
##
## Description:
## Calculate a series of simultaneous credible bounds from a samples matrix.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 02/10/2009   remove superfluous grid argument, add some checks,
##              really discard whole vectors and derive the SCB from the convex hull
##              of the rest.
#####################################################################################

## all methods assume that samples is a m by n matrix where
## m is the number of samples and n the number of parameters
## ==> hence each sample is a row in the matrix samples!

## scrHpd calculates a series of simultaneous credible bounds,
## minimizing the absolute distances to the mode vector at each gridpoint

scrHpd <- function(samples,          # sample matrix
                   mode = apply (samples, 2, median), # mode vector of length = ncol (samples)
                   level = 0.95     # credible level
                   )
{
    ## extracts
    nPars <- ncol(samples)                  # the number of parameters
    nSamples <- nrow(samples)                  # the number of samples

    ## checks
    if (nPars != length (mode))
        stop ("mode vector must have same length as samples matrix!")
    stopifnot(level > 0 && level < 1)

    ## absolute distance from mode vector
    distance <- abs (sweep (samples, 2, mode)) 

    ## Calculate a simultaneous (k/ nSamples)*100% credible band
    ## using the ranks approach:    
    k <- floor (level * nSamples)

    ## colwise (= elementwise) ranks of distances
    rankdistance <- apply (distance, 2, rank) 

    ## maximum ranks in each multivariate sample
    tstari <- apply (rankdistance, 1, max)

    ## sort the maximum ranks
    ordtstari <- sort.int (tstari, method = "quick") 

    ## the required rank, which divides the samples in contained and rejected samples for the SCB,
    ## is: 
    tstar <- ordtstari[k]               

    ## now which vectors are inside the SCB?
    whichInside <- tstari <= tstar
    ## note that sum(whichInside) is possibly larger than k, so we have
    ## a larger empirical coverage of the resulting SCB.

    ## reduce the samples matrix accordingly.
    samples <- samples[whichInside, ]

    ## the parameterwise ranges of these vectors form the SCB
    ## (just the convex hull of all sample vectors!) 
    ret <- apply(samples, 2, range)
    rownames(ret) <- c("lower", "upper")
    
    ## old code: it is not clear if here really a simultaenous credible band results!!
    ## because not whole sample vectors are discarded
    
    ## selectMat <- rankdistance <= tstar

    ## ret <- matrix (nrow = 2, ncol = nPars)
    ## for (i in seq_len (nPars)){
    ##     ret[, i] <- range (samples[selectMat[,i], i])
    ## }

    ## finally return the 2 x nPars matrix
    return (ret)
}
