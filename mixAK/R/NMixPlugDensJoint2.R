##
##  PURPOSE:   Computation of the pairwise joint densities
##             (plug-in version with supplied posterior summaries of mixture components)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTIONS: NMixPlugDensJoint2 (28/05/2009)
##
## ======================================================================

## *************************************************************
## NMixPlugDensJoint2
## *************************************************************
NMixPlugDensJoint2 <- function(x, ...)
{
  UseMethod("NMixPlugDensJoint2")
}  
