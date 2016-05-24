##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPredCDFMarg
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS: Y2T.NMixPredCDFMarg (09/06/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredCDFMarg
## *************************************************************
Y2T.NMixPredCDFMarg <- function(x, itrans=exp, ...)
{
  ### Untrans the grids and get T = itrans(Y)  
  if (length(itrans) == 1){  
    x$x <- lapply(x$x, itrans)
  }else{
    if (length(itrans) != length(x$x)) stop(paste("itrans must be of length ", length(x$x), sep=""))
    for (i in 1:length(x$dens)){
      x$x[[i]] <- itrans[[i]](x$x[[i]])
    }    
  }  
    
  return(x)
}  
