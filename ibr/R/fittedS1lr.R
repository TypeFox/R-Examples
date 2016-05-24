fittedS1lr <- function(n,U,tUy,eigenvaluesS1,ddlmini,k,rank) {
    if (rank<n) valpr0 <- eigenvaluesS1[-c(1:ddlmini,(rank+1):n)] else
        valpr0 <- eigenvaluesS1[-(1:ddlmini)]
    valpr <- rep(1,rank)
    valpr[-(1:ddlmini)] <- 1-(1-valpr0)^k
    fk <- U%*%(valpr*tUy[1:rank])
    return(list(fit=as.vector(fk),trace=sum(valpr)))
}
