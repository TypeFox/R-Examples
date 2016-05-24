#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getLogMargLik.R] by DSB Don 01/10/2009 17:35 (CEST)>
##
## Description:
## Extract log marginal likelihood from a model.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 03/09/2008   adapt for hyper-g methodology
## 04/09/2008   remove sst argument
## 19/10/2008   bugfix: R2 is inside the list, so we need double braces to get it.
#####################################################################################

getLogMargLik <- function (x,       # a valid BayesMfp-Object of length 1 (otherwise only first element recognized)
                           design=getDesignMatrix(x), # (centered) design matrix
                           nObs = nrow(design),
                           dim = ncol(design)
                           )
{
    ## select only the first element
    x <- x[1]

    ## todo: compute the R2 here! otherwise it is pretty unsafe...
    
    ## gather remaining arguments for logMargLik
    R2 <- x[[1]]$R2
    alpha <- attr(x, "priorSpecs")$a
    sst <- attr(x, "SST")
    
    ret <-
        .Call ("logMargLik", ## PACKAGE = "bfp",
               as.double(R2),
               as.integer(nObs),
               as.integer(dim),
               as.double(alpha),
               as.double(sst)
               )
    
    return (ret)
}
