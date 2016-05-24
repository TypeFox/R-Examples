##
##  PURPOSE:   Convert fitted distribution of Y=trans(T) into distribution of T=itrans(Y)
##             * method for objects of class NMixPredDensJoint2
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   17/08/2009
##
##  FUNCTIONS: Y2T.NMixPredDensJoint2 (17/08/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredDensJoint2
## *************************************************************
Y2T.NMixPredDensJoint2 <- function(x, itrans=exp, dtrans=function(x){return (1 / x)}, ...)
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
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      NAAM <- paste(i, "-", j, sep="")
      Jac <- matrix(rep(dx[[i]], length(dx[[j]])), nrow=length(dx[[i]]), ncol=length(dx[[j]])) * matrix(rep(dx[[j]], each=length(dx[[i]])), nrow=length(dx[[i]]), ncol=length(dx[[j]]))
      x$dens[[NAAM]] <- x$dens[[NAAM]] * Jac
      if (!is.null(x$densK)){
        for (k in 1:length(x$densK[[NAAM]])){
          x$densK[[NAAM]][[k]] <- x$densK[[NAAM]][[k]] * Jac
        }
      }          
    }  
  }  
  
  return(x)
}  
