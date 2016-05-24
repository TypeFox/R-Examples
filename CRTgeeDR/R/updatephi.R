#Exactly the function 'updatePhi' included in the package 'geeM',
#authored by Lee S. McDaniel and Nick Henderson
#under the GPL-2 license.
### Simple moment estimator of dispersion parameter
updatePhi <- function(Y, mu, VarFun, p, StdErr, included, includedlen){
  nn <- sum(includedlen)
  resid <- diag(StdErr %*% included %*% Diagonal(x=Y-mu))
  phi <- (1/(nn-p))*crossprod(resid, resid) 
  return(as.numeric(phi))	
}
