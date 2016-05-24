##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPlugCondDensMarg
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/08/2009
##
##  FUNCTIONS: Y2T.NMixPlugCondDensMarg (17/08/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPlugCondDensMarg
## *************************************************************
Y2T.NMixPlugCondDensMarg <- function(x, itrans=exp, dtrans=function(x){return (1 / x)}, ...)
{
  return(Y2T.NMixPredCondDensMarg(x=x, itrans=itrans, dtrans=dtrans, ...))
}
