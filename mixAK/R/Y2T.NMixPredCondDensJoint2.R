##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPredCondDensJoint2
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/08/2009
##
##  FUNCTIONS: Y2T.NMixPredCondDensJoint2 (17/08/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredCondDensJoint2
## *************************************************************
Y2T.NMixPredCondDensJoint2 <- function(x, itrans=exp, dtrans=function(x){return (1 / x)}, ...)
{
  p <- length(x$x)         ### number of margins (is 2 or more)

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
    for (i in 1:p){
      x$x[[i]] <- itrans[[i]](x$x[[i]])
      dx[[i]] <- dtrans[[i]](x$x[[i]])
    }  
  }      

  ### Multiply each bivariate density by appropriate jacobian (product of appropriate dx's)
  for (i in 1:length(x$dens)){           ## loop over values by which we condition
    for (j in 1:(p-1)){                  ## loop over margin 1
      if (j == x$icond) next             ## x$dens[[j-*]] is NA
      for (k in (j+1):p){                ## loop over margin 2
        if (k == x$icond) next           ## x$dens[[j-k]] is NA
        NAAM <- paste(j, "-", k, sep="")
        Jac <- matrix(rep(dx[[j]], length(dx[[k]])), nrow=length(dx[[j]]), ncol=length(dx[[k]])) * matrix(rep(dx[[k]], each=length(dx[[j]])), nrow=length(dx[[j]]), ncol=length(dx[[k]]))
        x$dens[[i]][[NAAM]] <- x$dens[[i]][[NAAM]] * Jac                
      }  
    }  
  }  
  
  return(x)
}
