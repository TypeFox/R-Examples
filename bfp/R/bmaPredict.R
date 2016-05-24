#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[bmaPredict.R] by DSB Don 17/06/2010 15:50 (CEST)>
##
## Description:
## BMA prediction for new data points.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 05/09/2008   correct: centering of x matrix, check ranges
##              numerical cancellations?! (in the summation)
## 10/11/2008   don't ask for the response when creating the model matrix
## 13/11/2008   use new internal function for creating model matrix for newdata
## 05/10/2009   some comments
## 17/06/2010   correct to use the shifts of the original data
#####################################################################################

## this is not a predict method for BmaSamples!
## it is a different predict "method" for BayesMfp

bmaPredict <- function (BayesMfpObject, # models over which to average the predictions
                        postProbs = posteriors (BayesMfpObject),
                                        # vector of posterior probabilites that will be normalized within
                        newdata        # new data as data.frame
                        )
{
    ## check that the probabilities are parallel to the models
    if (length (postProbs) != length (BayesMfpObject))
        stop ("postProbs has wrong length")

    ## construct the new uncentered data matrix
    tempX <- constructNewdataMatrix(BayesMfpObject=BayesMfpObject,
                                    newdata=newdata)

    ## and compute mean fit with respective original posterior coefficient mode for every model
    fitmat <- matrix (nrow = length (BayesMfpObject), ncol = nrow (tempX))
    for (i in seq_along (BayesMfpObject)){
        tempMod <- BayesMfpObject[i]

        ## get design matrix for the old data
        origDesign <- getDesignMatrix(tempMod, center=TRUE)

        ## and posterior parameter means
        post <- getPosteriorParms (BayesMfpObject[i])

        ## then change to the new uncentered covariates data matrix
        attr (tempMod, "x") <- tempX

        ## and use the old shifts to compute the
        ## correct design matrix for the new data
        tempDesign <- getDesignMatrix(tempMod, center=FALSE)
        tempDesign <- sweep(tempDesign,
                            MARGIN=2L,
                            attr(origDesign, "shifts"))

        ## to compute the fit
        fitmat[i,] <- tempDesign %*% post$mStar
    }

    ## average with probabilites as weights
    fitmat <- fitmat * (postProbs / sum (postProbs))
    return (colSums (fitmat))
}
