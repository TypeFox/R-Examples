#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getPostExpectedgAndShrinkage.R] by DSB Don 26/02/2009 12:31 (CET) on daniel@puc-home>
##
## Description:
## Extract/Compute posterior expected g and shrinkage factor from a model.
##
## History:
## 26/02/2009   file creation
#####################################################################################

getPostExpectedg <- function (x, # a valid BayesMfp-Object of length 1 (otherwise only first element
                                 # recognized) 
                              design=getDesignMatrix(x), # (centered) design matrix
                              nObs = nrow(design),
                              dim = ncol(design)
                              )
{
    ## select only the first element
    x <- x[1]

    ## gather remaining arguments for postExpectedg
    R2 <- x[[1]]$R2
    alpha <- attr(x, "priorSpecs")$a
    
    ret <-
        .Call ("postExpectedg", ## PACKAGE = "bfp",
               as.double(R2),
               as.integer(nObs),
               as.integer(dim),
               as.double(alpha)
               )
    
    return (ret)
}

getPostExpectedShrinkage <- function (x, # a valid BayesMfp-Object of length 1 (otherwise only
                                        # first element recognized) 
                                      design=getDesignMatrix(x), # (centered) design matrix
                                      nObs = nrow(design),
                                      dim = ncol(design)
                                      )
{
    ## select only the first element
    x <- x[1]

    ## gather remaining arguments for postExpectedShrinkage
    R2 <- x[[1]]$R2
    alpha <- attr(x, "priorSpecs")$a
    
    ret <-
        .Call ("postExpectedShrinkage", ## PACKAGE = "bfp",
               as.double(R2),
               as.integer(nObs),
               as.integer(dim),
               as.double(alpha)
               )
    
    return (ret)
}
