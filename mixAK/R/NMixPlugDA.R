##
##  PURPOSE:    Simple discriminant analysis based on mixture model
##              Determine component probabilities for "new" observations
##              based on posterior summary statistics of fitted components
##              (based on NMixMCMC)
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   16/06/2009  as NMixClust
##             22/11/2009  changed to NMixPlugDA
##
##  FUNCTION:  NMixPlugDA
##             
##
## ======================================================================

## *************************************************************
## NMixPlugDA
## *************************************************************
NMixPlugDA <- function(object, y)
{
  if (class(object) != "NMixMCMC") stop("object must be of class NMixMCMC")
  if (object$prior$priorK != "fixed") stop("number of mixture components was not fixed")

  if (object$nx_w > 1) stop("This function has not (yet) been implemented if a factor covariate on mixture weights is present.")
  
  if (object$dim == 1){
    if (is.data.frame(y) | is.matrix(y)){
      if (ncol(y) != 1) stop("y must have one column")
      namey <- rownames(y)
      laby <- colnames(y)[1]      
      y <- y[,1]
    }else{
      namey <- names(y)
      laby <- as.character(substitute(y))
      if (length(laby) > 1) laby <- "y"    ### this happens if, e.g., y=data[,"var"]
    }
    n <- length(y)
    RET <- data.frame(y)
    colnames(RET) <- laby
    
    W <- object$poster.mean.w
    Mu <- object$poster.mean.mu[,1] * object$scale$scale + object$scale$shift
    Var <- numeric(object$prior$Kmax)
    for (k in 1:object$prior$Kmax) Var[k] <- object$scale$scale^2 * as.numeric(object$poster.mean.Sigma[[k]])
    SD <- sqrt(Var)

    W <- matrix(rep(W, n), nrow=n, byrow=TRUE)
    Mu <- matrix(rep(Mu, n), nrow=n, byrow=TRUE)
    SD <- matrix(rep(SD, n), nrow=n, byrow=TRUE)    
    yy <- matrix(rep(y, object$prior$Kmax), nrow=n)
    ff <- dnorm(yy, Mu, SD)    
  }else{    
    if (!is.data.frame(y) & !is.matrix(y)) stop("y must be either a matrix or data.frame")
    if (is.matrix(y)){
      if (is.null(colnames(y))) colnames(y) <- paste("y", 1:object$dim, sep="")
      y <- as.data.frame(y)
    }
    namey <- rownames(y)    
    n <- nrow(y)
    RET <- y 
    
    W <- object$poster.mean.w
    ff <- matrix(NA, nrow=n, ncol=object$prior$Kmax)
    for (k in 1:object$prior$Kmax){
      Muk  <- object$poster.mean.mu[k,] * object$scale$scale + object$scale$shift
      Vark <- diag(object$scale$scale) %*% object$poster.mean.Sigma[[k]] %*% diag(object$scale$scale)
      ff[,k] <- dMVN(y, mean=Muk, Sigma=Vark)
    }     
  }

  Numer <- W*ff
  Denom <- apply(Numer, 1, sum)
  Poster <- Numer / matrix(rep(Denom, object$prior$Kmax), nrow=n)
  colnames(Poster) <- paste("prob", 1:object$prior$Kmax, sep="")
  rownames(Poster) <- namey

  ##RET <- cbind(RET, Poster)                           ### commented on 15/02/2010
  RET <- as.data.frame(Poster)                          ### added on 15/02/2010
  RET[,"component"] <- apply(Poster, 1, which.max)
  return(RET)
}
  
