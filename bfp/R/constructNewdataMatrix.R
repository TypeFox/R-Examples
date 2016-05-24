#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[constructNewdataMatrix.R] by DSB Mon 05/10/2009 10:01 (CEST)>
##
## Description:
## Internal function to construct the model matrix for new data based on the formula
## and scaling info in an existing BayesMfpObject, for passing it to prediction functions. 
##
## History:
## 13/11/2008   file creation
## 05/10/2009   comments
#####################################################################################


constructNewdataMatrix <- function(BayesMfpObject, # from BayesMfp
                                   newdata)        # new data as data.frame
{
    ## extract the covariates' formula and the names of the FP terms
    covariatesOnlyFormula <- attr (BayesMfpObject, "formula")[-2]
    bfpTermNames <- attr (BayesMfpObject, "termNames")$bfp

    ## get model matrix with the design variables for the new data
    newTerms <- terms (covariatesOnlyFormula, data = newdata) 
    ret <- model.matrix (newTerms, data = newdata)

    ## extract the transformation parameters
    tr <- attr (BayesMfpObject, "shiftScaleMax")

    ## and scale the new FP columns like for the old data
    for (bfp in bfpTermNames)           
    {
        ## transform
        ret[, bfp] <- (ret[, bfp] + tr[bfp, "shift"]) / tr[bfp, "scale"]

        ## and afterwards check range: all transformed numbers must be positive!
        if(any(ret[, bfp] <= 0))
            stop(simpleError(paste("New data for covariate", bfp, "is out of range.")))
    }

    ## return the correctly transformed model matrix
    return(ret)
}



