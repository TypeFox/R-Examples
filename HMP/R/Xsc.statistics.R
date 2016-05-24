Xsc.statistics <-
function(pi1, theta, nreads.data, pi0){
nReads <- sum(nreads.data)

Bj <- ((theta*(sum(nreads.data^2)-nReads)+nReads) / (nReads)^2) * (diag(pi0)-pi0 %*% t(pi0))
Xsc <- t(pi1-pi0) %*% MASS::ginv(Bj, tol=sqrt(.Machine$double.eps)) %*% (pi1-pi0)

return(Xsc)
}
