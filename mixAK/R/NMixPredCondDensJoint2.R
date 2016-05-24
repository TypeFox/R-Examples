##
##  PURPOSE:   Computation of all predictive pairwise bivariate conditional densities (given one margin)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   01/06/2009
##
##  FUNCTIONS: NMixPredCondDensJoint2 (01/06/2009)
##
## ======================================================================

## *************************************************************
## NMixPredCondDensJoint2
## *************************************************************
NMixPredCondDensJoint2 <- function(x, ...)
{
  UseMethod("NMixPredCondDensJoint2")
}  
