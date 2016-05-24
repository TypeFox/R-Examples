#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs with Hyper-g-prior
## Time-stamp: <[transformMfp.R] by DSB Fre 18/06/2010 14:04 (CEST)>
##
## Description:
## Transform a fitted mfp model into a BayesMfp model with the correct powers etc.
## to compare with other (true) BayesMfp models on the same (!) data.
##
## History:
## 26/02/2009   file creation: modify existing code from bfp/OzoneValidation.R
## 18/06/2010   sort fpNumeric accordingly to the order in BayesMfpObject, otherwise
##              the design matrix will not be correct!!
#####################################################################################

transformMfp <- function(mfpObject,
                         BayesMfpObject) # BayesMfp object, from which the first model is used for
                                        # imputation of the powers from mfpObject
{
    ## check distribution
    with(mfpObject$fit$family,
         stopifnot(identical(family, "gaussian"),
                   identical(link, "identity")))
    
    ## get the powers
    fpCharacter <- as.matrix(mfpObject$fptable)
    fpCharacter <- fpCharacter[, grep("power",
                                      colnames(fpCharacter))]

    ## transform to numbers
    oldOptions <- options(warn = -1)    # supresses warning for creation of NAs through as.numeric
    on.exit(options(oldOptions))        # transformation   
    
    fpNumeric <- matrix(data=as.numeric(fpCharacter),
                        nrow=nrow(fpCharacter),
                        ncol=ncol(fpCharacter),
                        dimnames=dimnames(fpCharacter))
    fpNumeric <- lapply(apply(fpNumeric,     
                              1,
                              sort,        # this sorts 
                              na.last=NA), # and discards NAs
                        as.vector)         # also get rid of attributes

    ## start return object
    ret <- BayesMfpObject[1]

    ## check equality of power names
    stopifnot(setequal(names(ret[[1]]$powers),
                       names(fpNumeric)))

    ## and then sort fpNumeric accordingly
    fpNumeric <- fpNumeric[names(ret[[1]]$powers)]
    
    ## impute into BayesMfp shell
    ret[[1]]$powers <- fpNumeric
    ret[[1]]$R2 <- with(mfpObject,
                        1 - deviance / null.deviance)
    ret[[1]]$logM <- getLogMargLik(ret)
    ret[[1]]$logP <- getLogPrior(ret)
    ret[[1]]$postExpectedg <- getPostExpectedg(ret)
    ret[[1]]$postExpectedShrinkage <- getPostExpectedShrinkage(ret)
    
    ## posterior values: only normalized estimate is available
    ret[[1]]$posterior <-  c(posterior=
                             exp(with(ret[[1]],
                                      logM + logP) - attr(BayesMfpObject, "logNormConst")),
                             frequency=NA)
    
    return(ret)    
}



