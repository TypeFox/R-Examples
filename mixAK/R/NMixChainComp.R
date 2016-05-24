##
##  PURPOSE:   Function to provide chains for mixture components
##             being shifted and scaled into the original scale and also
##             being relabeled if this is requested
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   15/07/2013
##
##  FUNCTION:  NMixChainComp
##             
##
## ======================================================================

## *************************************************************
## NMixChainComp
## *************************************************************
NMixChainComp <- function(x, relabel = TRUE, param)
{
  UseMethod("NMixChainComp") 
}
