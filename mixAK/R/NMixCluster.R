##
##  PURPOSE:   Clustering based on the MCMC output.
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:                  20/08/2014
##
##  FUNCTIONS:  NMixCluster
##
## ======================================================================

## *************************************************************
## NMixCluster
## *************************************************************
NMixCluster <- function(object, ...)
{
  UseMethod("NMixCluster") 
}

