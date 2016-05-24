#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[BayesMfpMethods.R] by DSB Fre 14/01/2011 16:28 (CET)>
##
## Description:
## Additional convenience methods for BayesMfp class objects.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 03/09/2008   modify for Hyper-g situation:
##              - as.data.frame method should extract posterior expectations and
##              R^2 as well
## 09/11/2008   add comment to predict method: this is a plug-in prediction
##              without interval generation.
## 10/11/2008   don't ask for the response inside the predict method
## 13/11/2008   use new internal function for creating model matrix for newdata
##              and scale the x matrix in tempMod
## 29/11/2008   update for new model prior option
## 17/06/2010   correct predict.BayesMfp to use the shifts of the original data
## 14/01/2011   remove spaces around pr$modelPrior
#####################################################################################

'[.BayesMfp' <- function (x, ...)       
{
    y <- NextMethod("[")
    mostattributes (y) <- attributes (x)
    names (y) <- names (x)[...]
    class(y) <- oldClass(x)
    y
}

####################################################################################################

print.BayesMfp <- function (x, ...)
{
    cat ("------------------------------ BayesMfp-Output ------------------------------\n")
    cat (length (x), "multivariable fractional polynomial model(s) of total (visited/cached)",
         attr (x, "numVisited"), "for following covariates:\n\n")
    cat ("fixed:                ", paste(attr (x, "termNames")$fixed, collapse = ", "), "\n")
    cat ("uncertain fixed form: ", paste(attr (x, "termNames")$uc, collapse = ", "), "\n")
    cat ("fractional polynomial:", paste(attr (x, "termNames")$bfp, collapse = ", "), "\n")
    pr <- attr (x, "priorSpecs")
    cat ("\nPrior parameter was a =", pr$a, "for the hyper-g prior,",
         "\nand a", pr$modelPrior, "prior on the model space was used.",
         "\n")
}

####################################################################################################

fitted.BayesMfp <- function (object, # valid BayesMfp object, only first element will be recognized
                             design = getDesignMatrix (object), # design matrix
                             post = getPosteriorParms (object, design = design), # posterior
                                        # parameters. if this default is used, then this is the
                                        # posterior expected mean, marginalized over g.
                             ...        # unused
                             )
{
    if (is.null (fit <- object[[1]]$fitted)){ # user may save them in here
        fit <- as.vector (design %*% post$mStar)
    }

    return (fit)
}

####################################################################################################

residuals.BayesMfp <- function (object, # valid BayesMfp object, only first element will be recognized
                                ...     # passed to fitted method
                                )
{
    fit <- fitted (object, ...)
    return (attr (object, "y") - fit)
}

####################################################################################################

## do a plug-in prediction at new data points, without credible intervals.
## Intervals can be obtained by first running BmaSamples and then extracting with the predict method
## of the BmaSamples class.
predict.BayesMfp <- function (object, # valid BayesMfp object, only first element will be recognized
                              newdata, # new covariate data with exactly the names (and preferably
                                       # ranges) as before
                              ...      # unused
                              )
{
    tempMod <- mod <- object[1]
    
    tempX <- constructNewdataMatrix(BayesMfpObject=object,
                                    newdata=newdata)   
    attr (tempMod, "x") <- tempX

    ## get design matrix for the old data
    origDesign <- getDesignMatrix(mod, center=TRUE)

    ## use its shifts to compute the correct design matrix for the new data
    design <- getDesignMatrix(tempMod, center=FALSE)
    design <- sweep(design,
                    MARGIN=2L,
                    attr(origDesign, "shifts"))
    
    ## and compute fit with original posterior coefficient mode   
    post <- getPosteriorParms (mod)
    fit <- as.vector (design %*% post$mStar)

    return (fit)
}

####################################################################################################

as.data.frame.BayesMfp <- function (x,  # valid BayesMfp object
                                    row.names = NULL, # optional rownames, default is to keep
                                        # listnames
                                    ...,        # unused
                                    freq = TRUE # should empirical frequencies be given?
                                    )
{
    ## posterior probabilites:
    ret <- data.frame (posterior = posteriors (x, ind = 1))
    row.names (ret) <- if (!is.null (row.names)) row.names else names (x)
    if (!is.null (attr (x, "chainlength")) && freq)
        ret$frequency <- posteriors (x, ind = 2)

    ## log marginal likelihood and log prior
    ret$logMargLik <- sapply (x, "[[", "logM")
    ret$logPrior <- sapply (x, "[[", "logP")

    ## posterior expected g and shrinkage g/(1 + g)
    ret$postExpectedg <- sapply(x, "[[", "postExpectedg")
    ret$postExpectedShrinkage <- sapply(x, "[[", "postExpectedShrinkage")

    ## coefficient of determination
    ret$R2 <- sapply(x, "[[", "R2")
    
    ## uncertain fixed form covariates:
    for (i in seq_along (ucNames <- attr (x, "termNames")$uc)){
        ret[[ucNames[i]]] <- sapply (x, function (one) i %in% one$ucTerms)
    }
    ## fractional polynomial:
    for (fpName in attr (x, "termNames")$bfp){
        ret[[fpName]] <- sapply (x, function (one) paste(one$powers[[fpName]], collapse = ", "))
    }

    ret
}
