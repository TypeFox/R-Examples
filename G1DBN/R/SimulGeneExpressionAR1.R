## __________________________________________________________
##
## FUNCTION SimulGeneExpressionAR1
##
## This function generates multivariate time series according 
## to the following first order Auto-Regressive process,
## X(t) = A X(t-1) + B + e(t)
## where var(e(t)) follows a zero-centered multivariate gaussian 
## distribution whose variance matrix S is diagonal.
##
## INPUT
##         - A : a matrix (p x p)
##         - B : a column vector (p x 1)
##         - X0 : a column vector (p x 1)
##         - SigmaEps : a column vector (p x 1)
##         - n : length of the time serie
## OUTPUT
##         - Xt : a matrix (n x p)
## __________________________________________________________
##

SimulGeneExpressionAR1<-function(A,B,X0,SigmaEps,n){
  
  ## initialization
  p <- dim(A)[1];
  Xt <- matrix(0,p,n)
  Xt[,1] <- X0
  
  ## Xt generation
  for (i in 2:n){
    Xt[,i] = A %*% Xt[,(i-1)] + B + rnorm(p,0,SigmaEps)}

  ## Setting Xt as a time series data set
  Xt <- ts(t(Xt), start=0, end=(n-1))

  return(Xt)
}
