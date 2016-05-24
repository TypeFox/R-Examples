fittedS1 <- function(n,U,tUy,eigenvaluesS1,ddlmini,k) { 
  valpr0 <- eigenvaluesS1[-(1:ddlmini)]
  valpr <- rep(1,n)
  valpr[-(1:ddlmini)] <- 1-(1-valpr0)^k
  fk <- U%*%(valpr*tUy)
  return(list(fit=as.vector(fk),trace=sum(valpr)))
}
