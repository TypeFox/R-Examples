##
##  PURPOSE:   Wishart distribution
##             * random numbers generation
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   12/11/2007
##
##  FUNCTION:  rWishart (12/11/2007)
##              -> renamed to rWISHART on 06/08/2013 (version 3.4-1)
##                 to avoid conflicts with rWishart from package stats
##
## ======================================================================

## *************************************************************
## rWishart
## *************************************************************
rWISHART <- function(n, df, S)
{
  thispackage <- "mixAK"

  if (is.null(dim(S))) wdim <- 1
  else{
    wdim <- nrow(S)
    if (ncol(S) != wdim) stop("S must be a square matrix")
  }
  lW <- (wdim*(wdim+1))/2
  
  if (df <= wdim - 1) stop(paste("df must be > ", wdim-1, sep=""))

  Si <- chol(S)
  Si <- chol2inv(Si)
  Sitri <- Si[lower.tri(Si, diag=TRUE)]

  ## Sample
  SAMPLE <- .C("rWishart_R", W       = double(n*lW),
                             dwork   = double(2*wdim*wdim),
                             int     = integer(1),
                             nu      = as.double(df),
                             invS    = as.double(Sitri),
                             dim     = as.integer(wdim),
                             npoints = as.integer(n),    
               PACKAGE = thispackage)

  Imat <- diag(wdim)
  rowsI <- row(Imat)[lower.tri(row(Imat), diag=TRUE)]
  colsI <- col(Imat)[lower.tri(col(Imat), diag=TRUE) ] 
  naam <- paste("(", rowsI, ".", colsI, ")", sep="")
  
  if (n == 1){
    if (wdim > 1){
      W <- diag(wdim)
      W[lower.tri(W, diag=TRUE)] <- SAMPLE$W
      W[upper.tri(W, diag=FALSE)] <- t(W)[upper.tri(t(W), diag=FALSE)]
      SAMPLE$W <- W
      rownames(SAMPLE$W) <- rownames(S)
      colnames(SAMPLE$W) <- colnames(S)      
    }  
  }
  else{
    if (wdim == 1){
      names(SAMPLE$W) <- 1:n
    }
    else{
      SAMPLE$W <- matrix(SAMPLE$W, byrow=TRUE, ncol=lW, nrow=n)         
      rownames(SAMPLE$W) <- 1:n
      colnames(SAMPLE$W) <- naam
    }  
  }  

  return(SAMPLE$W)
}


