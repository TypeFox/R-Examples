##
##  PURPOSE:   Moore-Penrose pseudoinverse of a symmetric squared matrix
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   24/01/2008
##
##  FUNCTIONS:  MatMPpinv
##
## ======================================================================

## *************************************************************
## MatMPpinv
## *************************************************************
MatMPpinv <- function(A)
{
  thispackage <- "mixAK"
  
  if (is.null(dim(A))) return(1/A)

  p <- nrow(A)
  if (ncol(A) != p) stop("A must be a squared matrix")

  LTA <- A[lower.tri(A, diag=TRUE)]
  UTA <- t(A)[lower.tri(A, diag=TRUE)]
  if (any(LTA != UTA)) stop("At this moment, only implemented for symmetric matrices.")

  RES <- .C("MPpinvSP", MPA=as.double(LTA), work=double(p + p*p + 3*p), err=as.integer(0), p=as.integer(p),
            PACKAGE=thispackage)

  MPA <- diag(p)
  MPA[lower.tri(MPA, diag=TRUE)] <- RES$MPA
  MPA[upper.tri(MPA, diag=FALSE)] <- t(MPA)[upper.tri(MPA, diag=FALSE)]  

  rownames(MPA) <- rownames(A)
  colnames(MPA) <- colnames(A)

  return(MPA)
}


