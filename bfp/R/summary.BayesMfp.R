#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[summary.BayesMfp.R] by DSB Fre 04/12/2009 17:06 (CET)>
##
## Description:
## Calculate and print the summary of a BayesMfp object.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 03/09/2008   add shrinkage factor argument
## 04/09/2008   remove sst arguments in call to getLogMargLik
## 02/10/2009   drop Bayes Factor computations, as they were wrong and are not important,
##              also be more careful about the "chainlength" attribute -- it is
##              only present when model sampling has been done.
#####################################################################################

`summary.BayesMfp` <-
    function (object,   # valid BayesMfp-object
              level = 0.95, # credible level for coefficients HPD intervals (only for 1 model)
              table = TRUE, # should a data.frame of the models be included?
              shrinkage=NULL,           # shrinkage factor used in the computations, defaults to the
                                        # posterior expected shrinkage factor
              ...       # arguments for as.data.frame
              )
{
    ## copy some attributes into return list
    objectAttributes <- attributes (object)
    ret <- objectAttributes[c ("call", "numVisited", "termNames", "shiftScaleMax",
                               "inclusionProbs")]
    if("chainlength" %in% names(objectAttributes))
        ret$chainlength <- objectAttributes$chainlength

    ## get the models table?
    if (table)
        ret$dataframe <- as.data.frame (object, ...)

    ## compute the local inclusion probabilities
    ret$localInclusionProbs <- inclusionProbs (object)

    ## if we have more than 1 model:
    if ((ret$nModels <- length (object)) > 1){
        ## extract the model probabilities
        ret$postProbs <- t (sapply (object, function (one) one[["posterior"]]))
        if (nrow (ret$postProbs) == 1){
            ret$postProbs <- t (ret$postProbs)
            colnames (ret$postProbs) <- "exact"
        } else {
            colnames (ret$postProbs) <- c ("normalized", "frequency")
        }
    } else {
        ## we have only one model, so summarize this one.
        ret$level <- level

        design <- getDesignMatrix (object)
        p <- ncol (design)

        shrinkage <-
            if(is.null(shrinkage))
                ## take the default for the shrinkage factor:
                object[[1]]$postExpectedShrinkage
            else
                ## check if it is valid:
                stopifnot(0 < shrinkage && shrinkage < 1)
        ret$shrinkage <- shrinkage
        
        post <- getPosteriorParms (object, design=design, shrinkage=shrinkage)

        summaryMat <- matrix (nrow = p, ncol = 3L,
                              dimnames =
                              list (rownames (post$mStar),
                                    c ("mode", "HPDlower", "HPDupper"))
                              )
        summaryMat <- as.data.frame (summaryMat)

        ## means
        summaryMat[,1] <- post$mStar

        ## HPD boundaries
        sStarSqrt <- with (post, sqrt (bStar / aStar * diag (VStar)))
        quant <- qt (p = 1 - level / 2, df = 2 * post$aStar)
        summaryMat[, 2] <- post$mStar - sStarSqrt * quant
        summaryMat[, 3] <- post$mStar + sStarSqrt * quant

        ## save
        ret$summaryMat <- summaryMat

        ## HPD for regression variance
        sigma2Mode <- with (post, bStar / (aStar + 1))
        rightLimit <- 2 * with (post, qinvGamma ((1 + level) / 2, aStar, bStar))

        EPS <- sqrt (.Machine$double.eps)

        outerdens = function(h, a, b){  # computes prob mass _out_ of the interval
            lowerInt <- c(EPS, sigma2Mode)
            upperInt <- c(sigma2Mode, rightLimit + EPS)

            schnitt.l = uniroot(function(x){dinvGamma(x, a, b) - h}, interval = lowerInt)$root
            schnitt.u = uniroot(function(x){dinvGamma(x, a, b) - h}, interval = upperInt)$root
            return(c(pinvGamma(schnitt.l, a, b) + pinvGamma(schnitt.u, a, b, lower.tail = FALSE),
                     schnitt.l, schnitt.u)
                   )
        }

        hInt <- with (post, c(dinvGamma (rightLimit, aStar, bStar), dinvGamma(sigma2Mode, aStar, bStar) - EPS))
        resulth = uniroot (function (h) {outerdens(h, post$aStar, post$bStar)[1] - (1 - level)}, interval = hInt)$root

        sigma2Sum <- summaryMat[1, 1:3]
        rownames (sigma2Sum) <- "sigma^2"
        sigma2Sum[1, 1] <- sigma2Mode
        sigma2Sum[1, 2:3] <- with (post, outerdens (resulth, aStar, bStar))[2:3]

        ## save
        ret$sigma2Sum <- sigma2Sum
    }
    class (ret) <- "summary.BayesMfp"
    return (ret)
}


####################################################################################################

print.summary.BayesMfp <- function (x, ...)
{
    cat ("------------------------------ BayesMfp-Output ------------------------------\n")
    if (x$nModels > 1){           # more than 1 model
        cat ("Original call:\n")
        print (x$call)
        cat ("\n")
        cat (x$nModels, " multivariable fractional polynomial model(s) of total (visited/cached)",
            x$numVisited,
             ifelse (!is.null (chl <- x$chainlength), paste("\n(during", chl, "sampling jumps)"), ""),
             " for following covariates:\n\n",
             sep = "")
        cat ("fixed:                ", paste(x$termNames$fixed, collapse = ", "), "\n")
        cat ("uncertain fixed form: ", paste(x$termNames$uc, collapse = ", "), "\n")
        cat ("fractional polynomial:\n\n")
        print (x$shiftScaleMax)
        cat ("\nDistribution of posterior probabilities:\n")
        print(summary (x$postProbs))
        cat ("\nOverall inclusion probabilities for covariates in question:\n")
        print (x$inclusionProbs)
        cat ("\nInclusion probabilities for covariates in question for this model selection:\n")
        print (x$localInclusionProbs)
        if (!is.null (tab <- x$dataframe)){
            cat ("\nOverview:\n")
            print (tab)
        }
    } else {                            # 1 model
        cat ("Model specification:\n")
        print (x$dataframe)

        cat("Chosen shrinkage factor:\n")
        print(signif(x$shrinkage, 3))
        
        cat ("\nPosterior summaries for single coefficients (",
             as.integer(x$level * 100), "% HPD intervals):\n", sep = "")

        print (x$summaryMat)

        cat ("Posterior mode and ", as.integer(x$level * 100),
             "% HPD interval for the regression variance:\n", sep = "")
        print (x$sigma2Sum)
    }
}
