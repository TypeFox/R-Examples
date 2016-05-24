##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS: Y2T (09/06/2009)
##
## ======================================================================

## *************************************************************
## Y2T
## *************************************************************
Y2T <- function(x, ...)
{
  UseMethod("Y2T")
}  
