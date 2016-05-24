##
##  PURPOSE:   Computation of the predictive marginal (univariate) cumulative distribution functions
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS: NMixPredCDFMarg (09/06/2009)
##
## ======================================================================

## *************************************************************
## NMixPredCDFMarg
## *************************************************************
NMixPredCDFMarg <- function(x, ...)
{
  UseMethod("NMixPredCDFMarg")
}  


