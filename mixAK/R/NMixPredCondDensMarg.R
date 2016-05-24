##
##  PURPOSE:   Computation of all predictive univariate conditional densities (given one margin)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   31/05/2009
##
##  FUNCTIONS: NMixPredCondDensMarg (31/05/2009)
##
## ======================================================================

## *************************************************************
## NMixPredCondDensMarg
## *************************************************************
NMixPredCondDensMarg <- function(x, ...)
{
  UseMethod("NMixPredCondDensMarg")
}  
