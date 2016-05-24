RefSigmaG <- function(sigma,Beta,X,y,delta,tol=0.0001,maxit=100,nitmon)
{
# Fixed point algorithm for scale
nit  <- 1
repeat{
  sigmao <- sigma
  sigma  <- ( RefAve2G(sigmao,Beta,X,y,delta)*sigmao^2/0.5 )^0.5; d <- sigma-sigmao
  if (nit==maxit | abs(d)<tol) break
  if(nitmon) cat(nit,sigma,"\n")
  nit <- nit+1}
list(sigma=sigma,nit=nit)}

