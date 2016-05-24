#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[inclusionProbs.R] by DSB Die 01/12/2009 10:58 (CET)>
##
## Description:
## Compute (model averaged) posterior inclusion probabilites for the uncertain
## variables (including FP variables) based on a BayesMfpObject.
##
## History:
## 04/07/2008   copy from thesis function collection.
#####################################################################################

`inclusionProbs` <-
    function (            # compute posterior inclusion probabilites based on BayesMfpObject
              BayesMfpObject,
              postProbs = posteriors (BayesMfpObject, ind = 1)
              )
{
    postProbs <- postProbs / sum (postProbs)
    inds <- attr (BayesMfpObject, "indices")
    termNames <- attr (BayesMfpObject, "termNames")

    nams <- unlist (termNames[c ("bfp", "uc")])
    ret <- numeric (length (nams))
    names (ret) <- nams

    i <- 0

    for (j in seq_along (inds$bfp)){
        i <- i + 1
        present <- sapply (BayesMfpObject, function (one) as.logical (length (one$powers[[j]])))
        ret[i] <- sum (present * postProbs)
    }

    for (j in seq_along (inds$ucList)){
        i <- i + 1
        present <- sapply (BayesMfpObject, function (one) any (j == one$ucTerms))
        ret[i] <- sum (present * postProbs)
    }

    return (ret)
}

