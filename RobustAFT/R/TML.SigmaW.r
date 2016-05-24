TML.SigmaW <- function(X,y,delta,sigma,sigma.t,mui,mui.t,wgt,const,cl,cu,tol,maxit,nitmon)
{
# Fixed point algorithm for scale
nit <- 1
repeat {
  sigmao <- sigma
  sigma  <- (TML.Ave2W(X,y,delta,sigmao,sigma.t, mui, mui.t,wgt, cl,cu,const)*sigmao^2/const)^0.5
  d <- sigma-sigmao
  if (nit==maxit | abs(d)< tol ) break
  if(nitmon) cat(nit,sigma,"\n")
  nit <- nit+1}
list(sigma=sigma,nit=nit)}

