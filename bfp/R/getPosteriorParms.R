#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getPosteriorParms.R] by DSB Fre 05/09/2008 13:41 (CEST) on daniel@puc.home>
##
## Description:
## Extract updated posterior parameters for the normal inverse gamma
## distribution from a model, given a shrinkage factor.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 03/09/2008   additional parameter "shrinkage" and new formulas inside, this is now only
##              conditional on a fixed shrinkage factor!
##              But the default guarantees correct posterior and predictive expectations.
## 04/09/2008   better computations, added OLS of non-intercept effects
##              and cholesky root of XtX to return list, don't compute the inverse of VStar
## 05/09/2008   select all rows of design matrix, correct betaOLS computation
#####################################################################################

`getPosteriorParms` <-
function (x,       # a valid BayesMfp-Object, only first list element will be recognized
          shrinkage=x[[1]]$postExpectedShrinkage,  # shrinkage factor
                                        # (defaults to the posterior expected g/(1+g) in this model)
          design = getDesignMatrix (x) # (centered) design matrix for the model
          )
{
    ret <- list ()

    n <- nrow (design)
    dim <- ncol(design)
    
    ## aStar
    ret$aStar <- (n - 1) / 2

    ## parts that are present in the null model too
    ret$VStar <- matrix(data=0, nrow=dim, ncol=dim)
    ret$VStar[1, 1] <- 1 / n

    ret$mStar <- matrix(data=0, nrow=dim, ncol=1)
    ret$mStar[1, 1] <- attr(x, "yMean")

    ## if this is not the null model
    if(dim > 1){
        ## then the rest is filled in now:
        nonInterceptColumns <- seq(from=2, to=dim)
        X <- design[, nonInterceptColumns]
        
        ## cholesky decomposition of XtX
        XtX <- crossprod (X)
        ret$XtXroot <- chol(XtX)
        
        ## fill in VStar
        ret$VStar[nonInterceptColumns, nonInterceptColumns] <-
            shrinkage * chol2inv(ret$XtXroot)

        ## compute betaOLS
        tmp <- forwardsolve(l=ret$XtXroot, x=crossprod(X, attr(x, "y")),
                            upper.tri=TRUE, transpose=TRUE) 
        ret$betaOLS <- backsolve(r=ret$XtXroot, x=tmp)

        ## fill in mStar
        ret$mStar[nonInterceptColumns, 1] <- shrinkage * ret$betaOLS
    }

    ## bStar
    ret$bStar <- (1 - shrinkage * x[[1]]$R2) * attr(x, "SST") / 2
    
    return (ret)
}

