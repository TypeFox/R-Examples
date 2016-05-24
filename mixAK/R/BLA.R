##
##  PURPOSE:   Best linear approximation 
##             (theoretical least squares)
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   19/11/2007
##
##  FUNCTIONS:  BLA
##
## ======================================================================

## *************************************************************
## BLA
## *************************************************************
BLA <- function(mean=c(0, 0),  Sigma=diag(2))
{
  thispackage <- "mixAK"
  
  p <- length(mean)
  LTp.1 <- ((p-1)*p)/2
  
  if (p <= 1) stop("mean must be of length >= 2")
  if (is.null(dim(Sigma))) stop("Sigma must be a matrix")
  if (p != nrow(Sigma) | p != ncol(Sigma)) stop("mean and Sigma are not consistent")
  Sigma <- Sigma[lower.tri(Sigma, diag=TRUE)]
  
  RES <- .C("BLA", beta   =double(p*p),
                   sigmaR2=double(p),
                   L      =double(LTp.1),            
                   err    =integer(1),
                   mu     =as.double(mean),
                   Sigma  =as.double(Sigma),
                   p      =as.integer(p),
            PACKAGE=thispackage)

  if (RES$err) stop("Sigma is singular")

  beta <- matrix(RES$beta, nrow=p, byrow=TRUE)
  rownames(beta) <- names(RES$sigmaR2) <- paste("x", 1:p, sep="")
  colnames(beta) <- paste("beta", 0:(p-1), sep="")  
  RET <- list(beta=beta,
              sigmaR2=RES$sigmaR2)

  return(RET)  
}

