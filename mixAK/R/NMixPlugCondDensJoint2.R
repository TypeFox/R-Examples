##
##  PURPOSE:   Computation of the pairwise conditional densities given one margin
##             (plug-in version with supplied posterior summaries of mixture components)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   29/05/2009
##
##  FUNCTIONS: NMixPlugCondDensJoint2 (29/05/2009)
##
## ======================================================================

## *************************************************************
## NMixPlugCondDensJoint2
## *************************************************************
NMixPlugCondDensJoint2 <- function(x, ...)
{
  UseMethod("NMixPlugCondDensJoint2")
}  
