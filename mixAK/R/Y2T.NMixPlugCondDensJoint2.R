##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPlugCondDensJoint2
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/08/2009
##
##  FUNCTIONS: Y2T.NMixPlugCondDensJoint2 (17/08/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPlugCondDensJoint2
## *************************************************************
Y2T.NMixPlugCondDensJoint2 <- function(x, itrans=exp, dtrans=function(x){return (1 / x)}, ...)
{
  return(Y2T.NMixPredCondDensJoint2(x=x, itrans=itrans, dtrans=dtrans, ...))
}
