##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPredCondDensMarg
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/08/2009
##
##  FUNCTIONS: Y2T.NMixPredCondDensMarg (17/08/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredCondDensMarg
## *************************************************************
Y2T.NMixPredCondDensMarg <- function(x, itrans=exp, dtrans=function(x){return (1 / x)}, ...)
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

  ### Multiply each conditional/marginal density by appropriate jacobian (t^{-1})
  for (i in 1:length(x$dens)){           ## loop over values by which we condition
    for (k in 1:length(x$dens[[i]])){    ## loop over margins
      if (k == x$icond){                 ## margin by which we condition
        x$dens[[i]][[k]] <- x$dens[[i]][[k]] * dx[[k]][i]       ## marginal density of margin by which we condition
      }else{
        x$dens[[i]][[k]] <- x$dens[[i]][[k]] * dx[[k]]          ## k-th conditional density
      }  
    }              
  }  

  ### Multiply each pointwise posterior quantile by appropriate jacobian (t^{-1})
  if (!is.null(x$prob)){
    for (j in 1:length(x$prob)){
      qnaam <- paste("q", x$prob[j]*100, "%", sep="")
      for (i in 1:length(x[[qnaam]])){           ## loop over values by which we condition
        for (k in 1:length(x[[qnaam]][[i]])){    ## loop over margins
          if (k == x$icond){                 ## margin by which we condition
            x[[qnaam]][[i]][[k]] <- x[[qnaam]][[i]][[k]] * dx[[k]][i]       ## marginal density of margin by which we condition
          }else{
            x[[qnaam]][[i]][[k]] <- x[[qnaam]][[i]][[k]] * dx[[k]]          ## k-th conditional density
          }  
        }              
      }        
    }  
  }  
  
  return(x)
}  
