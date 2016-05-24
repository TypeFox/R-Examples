##
##  PURPOSE:   Data manipulation for NMixMCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    14/02/2010 (by taking sub-code originally included in NMixMCMC function)
##
##  FUNCTIONS:  NMixMCMCdata
##
## ================================================================================================

## *************************************************************
## NMixMCMCdata
## *************************************************************
##
NMixMCMCdata <- function(y0, y1, censor)
{
  if (missing(y0)) stop("y0 (data) must be given")
  if (is.null(dim(y0))){
    y0 <- matrix(y0, ncol=1, nrow=length(y0))
    
    if (missing(censor)) censor <- rep(1, nrow(y0))
    if (!is.null(dim(censor))) stop("y0 and censor mismatch (dimensions)")
    censor <- matrix(censor, ncol=1, nrow=length(censor))
    
    if (missing(y1)){
      if (any(censor == 3)) stop("y1 must be given when there are some censor=3 values")      
      y1 <- rep(0, nrow(y0))
    }  
    if (!is.null(dim(y1))) stop("y0 and y1 mismatch (dimensions)")
    y1 <- matrix(y1, ncol=1, nrow=length(y1))
  }else{
    if (is.data.frame(y0)) y0 <- as.matrix(y0)
    
    if (missing(censor)) censor <- matrix(rep(1, nrow(y0)*ncol(y0)), nrow=nrow(y0), ncol=ncol(y0))
    if (is.null(dim(censor))) stop("y0 and censor mismatch (dimensions)")
    if (nrow(censor) != nrow(y0) | ncol(censor) != ncol(y0)) stop("y0 and censor mismatch (dimensions)")
    if (is.data.frame(censor)) censor <- as.matrix(censor)

    if (missing(y1)){
      if (any(censor == 3)) stop("y1 must be given when there are some censor=3 values")      
      y1 <- matrix(rep(0, nrow(y0)*ncol(y0)), nrow=nrow(y0), ncol=ncol(y0))
    }  
    if (is.null(dim(y1))) stop("y0 and y1 mismatch (dimensions)")
    if (nrow(y1) != nrow(y0) | ncol(y1) != ncol(y0)) stop("y0 and y1 mismatch (dimensions)")
    if (is.data.frame(y1)) y1 <- as.matrix(y1)    
  }  
  
  n <- nrow(y0)
  ###if (n < 2) stop("n (number of observations) must be at least 2")     ### commented on 15/02/2010
  p <- ncol(y0)
  LTp <- p * (p + 1)/2
  Imat <- diag(p)
  rowsI <- row(Imat)[lower.tri(row(Imat), diag=TRUE)]
  colsI <- col(Imat)[lower.tri(col(Imat), diag=TRUE)] 
  naamLTp <- paste(".", rowsI, ".", colsI, sep="")      
  
  if (any(is.na(y0))) stop("NA's in y0 are not allowed")
  
  if (missing(censor)) censor <- matrix(1, nrow=n, ncol=p) 
  if (is.null(dim(censor))) stop("y0 and censor mismatch (dimensions)") 
  if (nrow(censor) != n | ncol(censor) != p) stop("y0 and censor mismatch (nrow/ncol)")
  if (any(is.na(censor))) stop("NA's in censor are not allowed")
  if (any(!(as.matrix(censor) %in% c(0, 1, 2, 3)))) stop("censor values must be from {0, 1, 2, 3}")
  are.Censored <- any(censor != 1)
  are.Right    <- any(censor == 0)
  are.Exact    <- any(censor == 1)
  are.Left     <- any(censor == 2)
  are.Interval <- any(censor == 3)  
  
  is.Interval <- (censor == 3)

  if (missing(y1)){
    if (any(is.Interval)) stop("y1 must be given when there are some censor=3 values")
    y1 <- matrix(0, nrow=n, ncol=p)
  }
  if (!any(is.Interval)) y1 <- matrix(0, nrow=n, ncol=p)  
  if (is.null(dim(y1))) stop("y0 and y1 mismatch (dimensions)")
  if (nrow(y1) != n | ncol(y1) != p) stop("y0 and y1 mismatch (nrow/ncol)")  
  if (any(is.na(y1[is.Interval]))) stop("NA in the upper limit of the observed interval indicated")
  if (any(y0[is.Interval] >= y1[is.Interval])) stop("y0 and y1 mismatch (y0 >= y1 for some interval-censored observation)")
  y1[!is.Interval] <- 0

  RET <- list(y0           = y0,
              y1           = y1,
              censor       = censor,
              n            = n,
              p            = p,
              LTp          = LTp,
              naamLTp      = naamLTp,
              are.Censored = are.Censored,
              are.Right    = are.Right,
              are.Exact    = are.Exact,
              are.Left     = are.Left,
              are.Interval = are.Interval,
              is.Interval  = is.Interval)
  return(RET)
}  
