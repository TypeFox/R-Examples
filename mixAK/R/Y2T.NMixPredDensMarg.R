##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPredDensMarg
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS: Y2T.NMixPredDensMarg (09/06/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredDensMarg
## *************************************************************
Y2T.NMixPredDensMarg <- function(x, itrans=exp, dtrans=function(x){return (1 / x)}, ...)
{
  if (length(itrans) == 1){
    if (length(dtrans) != 1) stop("dtrans must be of length 1")
    
    ### Untrans the grids and get T = itrans(Y)
    x$x <- lapply(x$x, itrans)

    ### Multiply each marginal density by appropriate jacobian
    for (i in 1:length(x$dens)){
      x$dens[[i]] <- x$dens[[i]] * dtrans(x$x[[i]])
      if (!is.null(x$densK)){
        for (k in 1:length(x$densK[[i]])){
          x$densK[[i]][[k]] <- x$densK[[i]][[k]] * dtrans(x$x[[i]])
        }
      }  
    }    
  }else{
    if (length(itrans) != length(x$x)) stop(paste("itrans must be of length ", length(x$x), sep=""))
    if (length(dtrans) != length(x$x)) stop(paste("dtrans must be of length ", length(x$x), sep=""))

    ### Untrans the grids and get T = itrans(Y)
    ### and multiply each marginal density by appropriate jacobian
    for (i in 1:length(x$dens)){
      x$x[[i]] <- itrans[[i]](x$x[[i]])
      x$dens[[i]] <- x$dens[[i]] * dtrans[[i]](x$x[[i]])
      if (!is.null(x$densK)){      
        for (k in 1:length(x$densK[[i]])){
          x$densK[[i]][[k]] <- x$densK[[i]][[k]] * dtrans[[i]](x$x[[i]])
        }
      }  
    }  
  }      
  
  return(x)
}  
