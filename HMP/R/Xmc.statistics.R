Xmc.statistics <-
function(group.parameter, pi0){
K <- length(pi0)
numGroups <- length(group.parameter)
Xscg <- vector("numeric", numGroups)

for(i in 1:numGroups){
pi <- group.parameter[[i]][1:K]
theta <- group.parameter[[i]][K+1]
P <- length(group.parameter[[i]])
nreads.data <- group.parameter[[i]][(K+2):P]
nReads <- sum(nreads.data)

Bj <- ((theta*(sum(nreads.data^2)-nReads)+nReads) / nReads^2) * (diag(pi0)-pi0 %*% t(pi0))
Xscg[i] <- t(pi-pi0) %*% MASS::ginv(Bj, tol=sqrt(.Machine$double.eps)) %*% (pi-pi0)
}

return(sum(Xscg))
}
