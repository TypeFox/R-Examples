#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[summary.BmaSamples.R] by DSB Fre 02/10/2009 17:12 (CEST)>
##
## Description:
## Calculate and print the summary of a BmaSamples object.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 04/09/2008   add summary for shrinkage factor
#####################################################################################

summary.BmaSamples <-
    function (object,
              level = 0.95, # credible level for coefficients credible intervals
              hpd = TRUE, # should emprical hpd intervals be used (default) or simple quantile-based?
              ...       # unused
              )
{
    ## copy some elements into return list
    ret <- object[c ("sampleSize", "modelData", "modelFreqs")]
    ret$intervalType <- ifelse (hpd, "HPD", "equitailed")
    ret$level <- level

    ## construct coefficients' summary matrix
    summaryMat <- matrix (nrow = 0, ncol = 4,
                          dimnames = list (NULL, c ("median", "mean", "lower", "upper"))
                          )

    summarizeSamples <- function (mat){ # convenience function which summarizes a single matrix
        p <- ncol (mat)
        ret <- matrix (nrow = p, ncol = 4,
                       dimnames = list (colnames (mat), c ("median", "mean", "lower", "upper")))

        ret[, 1] <- apply (mat, 2, median)
        ret[, 2] <- colMeans (mat)

        ret[, 3:4] <-
            if (hpd){
                t (apply (mat, 2, empiricalHpd, level = level))
            } else {
                alpha <- 1 - level
                t (apply (mat, 2, quantile, probs = c (alpha / 2, 1 - alpha / 2)))
            }

        return (ret)
    }

    ## process samples
    summaryMat <- rbind (summaryMat, summarizeSamples (object$fixed)) # fixed coefs
    for (ucMat in object$uc){
        summaryMat <- rbind (summaryMat, summarizeSamples (ucMat)) # uncertain covariates coefs
    }
    ## and save
    ret$summaryMat <- summaryMat

    ## summary for regression variance
    sigma2Sum <- summaryMat[1, 1:4, drop = FALSE]
    rownames (sigma2Sum) <- "sigma^2"
    sigma2Sum[1:2] <- with (object, c (median (sigma2), mean (sigma2)))
    sigma2Sum[3:4] <-
        if (hpd){
            empiricalHpd (object$sigma2, level = level)
        } else {
            alpha <- 1 - level
            quantile (object$sigma2, probs = c (alpha / 2, 1 - alpha / 2))
        }
    ## save
    ret$sigma2Sum <- sigma2Sum

    ## summary for shrinkage factor
    shrinkageSum <- summaryMat[1, 1:4, drop = FALSE]
    rownames (shrinkageSum) <- "shrinkage factor"
    shrinkageSum[1:2] <- with (object, c (median (shrinkage), mean (shrinkage)))
    shrinkageSum[3:4] <-
        if (hpd){
            empiricalHpd (object$shrinkage, level = level)
        } else {
            alpha <- 1 - level
            quantile (object$shrinkage, probs = c (alpha / 2, 1 - alpha / 2))
        }
    ## save
    ret$shrinkageSum <- shrinkageSum
    
    class (ret) <- "summary.BmaSamples"
    return (ret)

}

####################################################################################################

print.summary.BmaSamples <- function (x,
                                      table = TRUE, # should the model table been shown?
                                      ... # unused
                                      )
{
    cat ("----------------------------- BmaSamples-Output -----------------------------\n")
    cat (x$sampleSize, "samples of the Bayesian model average over",
         nrow (x$modelData), "BayesMfp models\n")

    if (table){
        cat ("\nModel overview (sample frequencies attached right):\n")
        print (x$modelData)
    }

    cat ("\nPosterior summaries for single coefficients (",
         as.integer(x$level * 100), "% ", x$intervalType, " intervals):\n", sep = "")
    print (x$summaryMat)

    cat ("\nPosterior mode and ", as.integer(x$level * 100),
         "% ", x$intervalType, " interval for the regression variance:\n", sep = "")
    print (x$sigma2Sum)

    cat ("\nPosterior mode and ", as.integer(x$level * 100),
         "% ", x$intervalType, " interval for the shrinkage factor:\n", sep = "")
    print (x$shrinkageSum)

}
