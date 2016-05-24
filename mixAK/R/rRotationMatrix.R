##
##  PURPOSE:   Generate random rotation matrices
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/01/2008
##
##  FUNCTIONS:  rRotationMatrix
##
## ======================================================================

## *************************************************************
## rRotationMatrix
## *************************************************************
rRotationMatrix <- function(n, dim)
{
  thispackage <- "mixAK"

  if (n <= 0) stop("n must be positive")
  if (dim <= 0) stop("dim must be positive")
  dim2 <- dim^2
  
  ldwork <- dim2 + dim + 2*dim + dim2

  PC <- .C("RotationMatrix_R",
           P=double(n*dim2), dwork=double(ldwork), pivot=integer(dim), err=as.integer(0), dim=as.integer(dim), n=as.integer(n),
           PACKAGE=thispackage)
  if (PC$err) warning("There might be problems")
  
  if (n == 1){
    return(matrix(PC$P, nrow=dim, ncol=dim))
  }else{
    P <- list()
    for (i in 1:n){
      P[[i]] <- matrix(PC$P[((i-1)*dim2+1):(i*dim2)], nrow=dim, ncol=dim)
    }
    return(P)
  }   
}  
