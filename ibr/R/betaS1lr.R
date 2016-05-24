betaS1lr <- function(n,U,tUy,eigenvaluesS1,ddlmini,k,lambda,rank,Rm1U,index0) {
  prov <- rev(sumvalpr(k,n,rev(1-eigenvaluesS1),n-index0+1,n-ddlmini+1))
  beta <- U%*%(prov[1:rank]*tUy[1:rank])
  beta <- t(eigenvaluesS1[1:rank]*t(Rm1U)) %*%t(U) %*% beta
  return(beta)
}
