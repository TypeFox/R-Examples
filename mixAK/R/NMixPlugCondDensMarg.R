##
##  PURPOSE:   Computation of all univariate conditional densities (given one margin)
##             (plug-in version with supplied posterior summaries of mixture components)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTIONS: NMixPlugCondDensMarg (28/05/2009)
##
## ======================================================================

## *************************************************************
## NMixPlugCondDensMarg
## *************************************************************
NMixPlugCondDensMarg <- function(x, ...)
{
  UseMethod("NMixPlugCondDensMarg")
}  
