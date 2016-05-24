##
##  PURPOSE:   Computation of the marginal densities
##             (plug-in version with supplied posterior summaries of mixture components)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTIONS: NMixPlugDensMarg (28/05/2009)
##
## ======================================================================

## *************************************************************
## NMixPlugDensMarg
## *************************************************************
NMixPlugDensMarg <- function(x, ...)
{
  UseMethod("NMixPlugDensMarg")
}  
