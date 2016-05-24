#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[BmaSamplesMethods.R] by DSB Fre 14/01/2011 16:28 (CET)>
##
## Description:
## Additional convenience methods for BmaSamples class objects.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 04/09/2008   modify for new hyper-g methodology
## 09/11/2008   add predict method to extract predictions for new data from
##              BmaSamples object
## 10/11/2008   and a corresponding print method
## 29/11/2008   update for new model prior option
## 14/01/2011   remove spaces around pr$modelPrior
#####################################################################################


print.BmaSamples <- function (x,
                              ...       # unused
                              )
{
    cat ("----------------------------- BmaSamples-Output -----------------------------\n")
    cat (x$sampleSize, "samples of the Bayesian model average over",
         nrow (x$modelData), "BayesMfp models",
         "\nwhich (partly) include the following covariates:\n\n")
    tn <- x$termNames
    cat ("fixed:                ", paste (tn$fixed, collapse = ", "), "\n")
    cat ("uncertain fixed form: ", paste (tn$uc, collapse = ", "), "\n")
    cat ("fractional polynomial:", paste (tn$bfp, collapse = ", "), "\n")
    pr <- x$priorSpecs
    cat ("\nPrior parameter was a =", pr$a, "for the hyper-g prior,",
         "\nand a", pr$modelPrior, "prior on the model space was used.",
         "\n")
}

####################################################################################################

fitted.BmaSamples <- function (object, ...)
{
    probs <- object$modelData$bmaProb
    weightedFits <- object$fitted * probs
    return (colSums (weightedFits))
}

####################################################################################################

residuals.BmaSamples <- function (object, ...)
{
    fit <- fitted (object, ...)
    return (object$y - fit)
}

####################################################################################################

## extract predictions and their intervals from a BmaSamples object
predict.BmaSamples <- function (object, # valid BmaSamples object
                                level=0.95, # credible level for predictions credible intervals
                                hpd = TRUE, # should emprical hpd intervals be used (default) or simple quantile-based?
                                ...         # unused
                                )
{
    ## copy some elements into return list
    ret <- list()

    ret$intervalType <- ifelse (hpd, "HPD", "equitailed")
    ret$level <- level
    ret$newdata <- object$newdata

    ret$sampleSize <- ncol(object$predictions)
    ret$nModels <- nrow(object$modelData)
    
    ## construct predictions summary matrix
    ret$summaryMat <- matrix (nrow = nrow(ret$newdata),
                              ncol = 4,
                              dimnames =
                              list (rownames(ret$newdata),
                                    c ("median", "mean", "lower", "upper")))                          
    
    ret$summaryMat[, "median"] <- apply (object$predictions,
                                         1,
                                         median)
    ret$summaryMat[, "mean"] <- rowMeans (object$predictions)

    ret$summaryMat[, c("lower", "upper")] <-
        if (hpd)
        {
            t (apply (object$predictions,
                      1,
                      empiricalHpd,
                      level = level))
        }
        else
        {
            alpha <- 1 - level
            t (apply (object$predictions,
                      1,
                      quantile,
                      probs =
                      c (alpha / 2,
                         1 - alpha / 2)))
        }
    
    class(ret) <- "predict.BmaSamples"
    return (ret)
}

####################################################################################################

print.predict.BmaSamples <- function (x,
                                      ... # unused
                                      )
{
    cat ("Prediction for new data points using Bayesian model averaging\n")
    cat ("over", x$nModels, "BayesMfp models, from which", x$sampleSize, "posterior predictive\n")
    cat("samples were drawn.")

    cat ("\nSummaries of these samples (",
         as.integer(x$level * 100), "% ", x$intervalType, " intervals):\n", sep = "")
    print (x$summaryMat)
}

