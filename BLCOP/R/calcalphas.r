###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# CAPMList
# Author: Francisco
###############################################################################
# DESCRIPTION: Computes the CAPM alphas of a set of returns with a market index and a give risk-free rate
# KEYWORDS: math
###############################################################################

CAPMList <- function
(
    returns,                            # matrix of returns
    marketIndex,                        # vector or time series of market index
    riskFree = NULL,                    # risk-free rate of return
    regFunc = BLCOPOptions("regFunc"),  # function to use to perform regression
    coeffExtractFunc = NULL,            # function to extract intercept and betas of regression
    ...                                 # additional parameters to the regression function
)
{
    
    CAPMfits <- vector(mode = "list", length = ncol(returns))
    regFunc <- match.fun(regFunc)
    # if risk-free rate is missing, replace it with 0
    if(is.null(riskFree))
        riskFree <- rep(0, length(marketIndex) )
    for(i in 1:ncol(returns))
        CAPMfits[[i]] <- regFunc(returns[,i] - riskFree ~ I(marketIndex - riskFree), ...)
    
    coeffs <- lapply(CAPMfits, coef)
    
    # Note: this is not a generic function as it does not work for all results
    # from linear regression functions
    if(is.null(coeffExtractFunc))    
    {
        coeffExtractFunc <- function(fit)
        {
            c(fit["(Intercept)"], fit["I(marketIndex - riskFree)"])
        }
    }
    results <- try(sapply(coeffs, coeffExtractFunc))
    if(inherits(results, "try-error"))
        stop("Unable to extract coefficients from regression results")
    else
    {  
      data.frame("alphas" = results["(Intercept)",], "betas" = results["I(marketIndex - riskFree)",], row.names = colnames(returns) )
        
    }
}     