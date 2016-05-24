##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPredCondCDFMarg
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   08/05/2010
##
##  FUNCTIONS: Y2T.NMixPredCondDensMarg (08/05/2010)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredCondCDFMarg
## *************************************************************
Y2T.NMixPredCondCDFMarg <- function(x, itrans=exp, dtrans=function(x){return (1 / x)}, ...)
{
  if (length(itrans) == 1){
    if (length(dtrans) != 1) stop("dtrans must be of length 1")
    
    ### Untrans the grids and get T = itrans(Y)
    x$x <- lapply(x$x, itrans)

    ### Compute derivatives in transformed grid
    dx <- lapply(x$x, dtrans)
  }else{
    if (length(itrans) != length(x$x)) stop(paste("itrans must be of length ", length(x$x), sep=""))
    if (length(dtrans) != length(x$x)) stop(paste("dtrans must be of length ", length(x$x), sep=""))

    ### Untrans the grids and get T = itrans(Y)
    ### Compute derivatives in transformed grid
    dx <- list()
    for (k in 1:length(x$x)){
      x$x[[k]] <- itrans[[k]](x$x[[k]])
      dx[[k]] <- dtrans[[k]](x$x[[k]])      
    }  
  }      
  
  return(x)
}  
