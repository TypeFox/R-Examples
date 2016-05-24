##
##  PURPOSE:   Computation of the predictive marginal (univariate) densities
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   03/12/2007
##
##  FUNCTIONS: NMixPredDensMarg (03/12/2007)
##
## ======================================================================

## *************************************************************
## NMixPredDensMarg
## *************************************************************
NMixPredDensMarg <- function(x, ...)
{
  UseMethod("NMixPredDensMarg")
}  


