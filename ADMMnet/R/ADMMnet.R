

####################################################
#####  ADMMnet (ADMM-L0)                       #####
#####  Penalty: L1, L2, Laplacian              #####
#####  Algorithm: one-step coordinate descent  #####
####################################################

ADMMnet=function(x, y, family=c("gaussian", "cox"), penalty=c("Lasso","Enet", "Net"), Omega=NULL, alpha=1.0, lambda=NULL, nlambda=50, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, adaptive=c(FALSE, TRUE), aini=NULL, isd=FALSE, keep.beta=FALSE, ifast=TRUE, thresh=1e-7, maxit=1e+5) {
  
  #fcall=match.call()
  family=match.arg(family)
  penalty=match.arg(penalty)
  if (penalty=="Lasso") {
    penalty="Enet"
    alpha=1.0
  }
  
  if (penalty=="Net" & is.null(Omega)) {
    penalty="Enet"
    cat("Enet was performed as no input of Omega")
  }
  
  if (family == "gaussian") {
    fit=switch(penalty,
               "Enet"=EnetLm(x,y,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,adaptive[1],aini,isd,keep.beta,thresh,maxit),
               "Net"=NetLm(x,y,Omega,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,adaptive,aini,isd,keep.beta,thresh,maxit))
    fit$family="gaussian"
    
  } else if (family == "cox") {
    fit=switch(penalty,
               "Enet"=EnetCox(x,y,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,adaptive[1],aini,isd,keep.beta,ifast,thresh,maxit),
               "Net"=NetCox(x,y,Omega,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,adaptive,aini,isd,keep.beta,ifast,thresh,maxit))
    fit$family="cox"
  }
  
  #fit$call=fcall
  class(fit)="ADMMnet"
  return(fit)  
}




