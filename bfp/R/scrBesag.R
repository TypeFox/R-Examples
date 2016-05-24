#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[scrBesag.R] by DSB Mit 26/01/2011 14:54 (CET)>
##
## Description:
## Calculate a series of simultaneous credible bounds from a samples matrix using
## the original Besag algorithm.
##
## History:
## 26/01/2011   copy from hypergsplines package, slightly modified to have
##              2 rbinded rows as result (instead of 2 cbinded columns) and
##              samples x parameters layout.
#####################################################################################

scrBesag <- function(samples,
                     level=0.95)
{
    n <- nrow(samples)
    p <- ncol(samples)

    samples.sorted <- matrix(0, n, p)
    samples.rank <- matrix(0, n, p)
    for(i in 1:p)
    {
        samples.sorted[,i] <- sort(samples[,i])
        samples.rank[,i] <- rank(samples[,i])
    }
    
    k <- trunc(0.95*n)+1
    helpset <- rep(0, n)
    for(i in 1:n)
        helpset[i] <- max(c(n+1-min(samples.rank[i,]), max(samples.rank[i,])))
    helpset <- sort(helpset)
    t.star <- helpset[k]

    lower <- samples.sorted[n+1-t.star,]
    upper <- samples.sorted[t.star,]

    return(rbind(lower=lower, upper=upper))
}
