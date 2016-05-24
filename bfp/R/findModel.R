#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[findModel.R] by DSB Mon 05/10/2009 10:24 (CEST)>
##
## Description:
## Find a model with specific powers and nontransformed uncertainty terms
## among the models in a BayesMfp oject.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 12/02/2009   faster implementation using "match".
#####################################################################################

`findModel` <-
    function (                              # returns the first index of
              model,                        # model in
              BayesMfpObject                # this model collection
              )                             # model must be a list with elements powers and ucTerms
                                        # at best, copy a model object and change the elements appropriately!
{
    ## the subset which is compared
    subset <- c("powers", "ucTerms")

    ## make sure the mode of the ucTerms is correct
    model$ucTerms <- as.integer(model$ucTerms)
    x <- list(model[subset])

    ## table with the model specifications
    table <- lapply(BayesMfpObject, "[", subset)

    ## the lookup
    match(x, table)
}

