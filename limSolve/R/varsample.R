
##==============================================================================
## varsample  : Samples variable equations
##==============================================================================

varsample <- function(X,  EqA, EqB=NULL)  {

  if (! is.matrix(EqA) & ! is.null(EqA))
    EqA <- t(as.matrix(EqA))
  if (! is.matrix(X) & ! is.null(X))
    X <- t(as.matrix(X))

  if (is.null(EqB))
    EqB <- rep(0,nrow(EqA))
  Var <- NULL
  if (ncol(X) != ncol(EqA))
    stop("matrix X and EqA not compatible")
  for (i in 1:nrow(X))
    Var <- rbind(Var,t(EqA%*%X[i,]-EqB))
  return(Var)
}
