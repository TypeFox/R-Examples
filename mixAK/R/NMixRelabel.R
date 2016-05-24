##
##  PURPOSE:   Re-label mixture components in the MCMC output obtained using the NMixMCMC/GLMM_MCMC function
##             The function also re-calculates 'poster.comp.prob_u' and 'poster.comp.prob_b'.
##             It is only applicable for models with fixed number of mixture components.
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:                  08/02/2010
##
##  FUNCTIONS:  NMixRelabel
##
## ======================================================================

## *************************************************************
## NMixRelabel
## *************************************************************
NMixRelabel <- function(object, type=c("mean", "weight", "stephens"), par, ...)
{
  UseMethod("NMixRelabel") 
}

