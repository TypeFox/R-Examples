#*** rMVNorm.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:  06/12/2006
##
## PURPOSE:  - random numbers from multivariate normal distribution
##           
## FUNCTIONS:   rMVNorm
##
#* ********************************************************************

##########################################################################################
### rMVNorm: Random number generation from the multivariate normal distribution
###
##########################################################################################
rMVNorm <- function(n, mean=0, Sigma=1, Q, param=c("standard", "canonical"))
{
  thispackage <- "bayesSurv"
  #thispackage <- NULL
  
  param <- match.arg(param)

  if (missing(Q)){
    if (param == "canonical") stop("Q must be given if param = canonical")
    if (missing(Sigma)) stop("Sigma must be given")
    QS <- Sigma
    version <- 0
  }
  else{
    QS <- Q
    if (param == "standard")  version <- 1
    if (param == "canonical") version <- 2
  }
    
  if (is.null(dim(QS))){
    if (length(QS) > 1) stop("Q and/or Sigma must be either a number of a square matrix")
    nx <- 1
  }  
  else{
    if (ncol(QS) != nrow(QS)) stop("Q and/or Sigma must be a square matrix")
    nx <- ncol(QS)
  }
  QStri <- QS[lower.tri(QS, diag=TRUE)]

  if (length(mean) == 1){
    mean <- rep(mean, nx)
    names(mean) <- paste("x", 1:nx, sep="")
  }
  else{
    if (length(mean) != nx) stop(paste("mean must be of length ", nx, sep=""))
    if (is.null(names(mean))) names(mean) <- paste("x", 1:nx, sep="")
  }
  rownames(QS) <- names(mean)

  ## Sample
  SAMPLE <- .C("rmvnormR2006", x=double(n*nx),
                               mu = as.double(mean),
                               QS = as.double(QStri),
                               err = integer(1),
                               nx = as.integer(nx),
                               nrandom = as.integer(n),
                               version = as.integer(version),
               PACKAGE = thispackage)
               
  if (n == 1){
    names(SAMPLE$x) <-  names(mean)
  }
  else{
    SAMPLE$x <- matrix(SAMPLE$x, byrow=TRUE, ncol=nx, nrow=n)
    colnames(SAMPLE$x) <- names(mean)
    rownames(SAMPLE$x) <- 1:n    
  }

  return(SAMPLE$x)
}


