fittedA <-  function(n,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,k){
  ## premier test
  prov <- rep(1,n)
  if (ddlmini>=1) {
    valpr0 <- 1-eigenvaluesA[-(1:ddlmini)]
    prov[-(1:ddlmini)] <- 1-valpr0^k
  } else {
    valpr0 <- 1-eigenvaluesA
    prov <- 1-valpr0^k
  }
  trace <- sum(prov)
  prov1 <- matrix(prov*as.vector(tPADmdemiY),n,1)
  return(list(fit=as.vector(DdemiPA%*%prov1),trace=trace))
}
