##
##  PURPOSE:   Computation of all predictive univariate conditional cdf's (given one margin)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   06/05/2010
##
##  FUNCTIONS: NMixPredCondCDFMarg (06/05/2010)
##
## ======================================================================

## *************************************************************
## NMixPredCondCDFMarg
## *************************************************************
NMixPredCondCDFMarg <- function(x, ...)
{
  UseMethod("NMixPredCondCDFMarg")
}  
