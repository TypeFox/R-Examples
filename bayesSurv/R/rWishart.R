#*** rWishart.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:  11/12/2006
##
## PURPOSE:  - random numbers from the Wishart distribution
##           
## FUNCTIONS:   rWishart
##
#* ********************************************************************

##########################################################################################
### rWishart: Random number generation from the Wishart distribution
###
##########################################################################################
rWishart <- function(n, df, S)
{
  thispackage <- "bayesSurv"
  #thispackage <- NULL

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
  SAMPLE <- .C("rwishartR3", W    = double(n*lW),
                             work = double(2*wdim*wdim),
                             df   = as.double(df),
                             invS = as.double(Sitri),
                             dim  = as.integer(wdim),
                             nrandom = as.integer(n),    
               PACKAGE = thispackage)

  Imat <- diag(wdim)
  rowsI <- row(Imat)[lower.tri(row(Imat), diag=TRUE)]
  colsI <- col(Imat)[lower.tri(col(Imat), diag=TRUE) ] 
  naam <- paste("(", rowsI, ".", colsI, ")", sep="")
  
  if (n == 1){
    if (wdim > 1){
      names(SAMPLE$W) <- naam
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

