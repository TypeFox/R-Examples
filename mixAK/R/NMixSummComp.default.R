##
##  PURPOSE:   MCMC for a normal mixture model
##             * basic posterior summary for mixture components (after re-labeling to achieve
##               at least some identifiability) in a model with fixed number of components
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/06/2010
##
##  FUNCTION:  NMixSummComp.default
##             
##
## ======================================================================

## *************************************************************
## NMixSummComp
## *************************************************************
NMixSummComp.default <- function(x)
{
  stop("There is no default NMixSummComp method.")  
}
