##
##  PURPOSE:   Wishart distribution
##             * (log-)density
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   12/11/2007
##
##  FUNCTIONS:  dWishart (16/01/2008)
##              -> renamed to dWISHART on 06/08/2013 (version 3.4-1)
##                 to avoid conflicts with rWishart from package stats
##              
##
## ======================================================================

## *************************************************************
## dWishart
## *************************************************************
dWISHART <- function(W, df, S, log=FALSE)
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

  if (is.null(dim(W))){
    if (wdim == 1) n <- length(W)
    else           stop("W must be a matrix")
  }else{
    if (nrow(W) == wdim & ncol(W) == wdim){
      n <- 1
      W <- W[lower.tri(W, diag=TRUE)]
    }else{
      if (ncol(W) != lW) stop(paste("W must have ", lW, " columns (lower triangles of sampled W form rows of the argument W)", sep=""))
      n <- nrow(W)
      W <- as.numeric(t(W))
    }  
  }  

  ## Compute log-density
  lDens <- .C("ldWishart_R", ldens            = double(n),
                             W.L              = double(n*lW),
                             log.sqrt.detW    = double(n), 
                             log.const        = double(1),
                             invS.L           = double(lW),
                             log.sqrt.detinvS = double(1),
                             err              = as.integer(0),
                             W                = as.double(W),
                             nu               = as.double(df),
                             invS             = as.double(Sitri),
                             dim              = as.integer(wdim),
                             npoints          = as.integer(n),              
              PACKAGE = thispackage)
  
  if (!log) lDens$ldens <- exp(lDens$ldens)
  return(lDens$ldens)  
}


