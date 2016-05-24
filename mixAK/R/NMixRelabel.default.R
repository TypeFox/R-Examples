##
##  PURPOSE:   Re-labeling of the MCMC output.
##             * default method
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   08/02/2010
##
##  FUNCTION:  NMixRelabel.default (08/02/2010) 
##
## ======================================================================

## *************************************************************
## NMixRelabel.default
## *************************************************************
NMixRelabel.default <- function(object, type=c("mean", "weight", "stephens"), par, ...)
{
  stop("There is no default NMixRelabel method.")
}  
