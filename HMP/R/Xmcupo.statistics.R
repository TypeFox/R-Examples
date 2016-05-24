Xmcupo.statistics <-
function(group.parameter, K){
numGrps <- length(group.parameter)
Xscg <- matrix(0, K+1, numGrps)
for(i in 1:numGrps){
theta <- group.parameter[[i]][K+1]
nreads.data <- group.parameter[[i]][-(1:(K+1))]
nReads <- sum(nreads.data)

N_1Cj <- (theta*(sum(nreads.data^2)-nReads) + nReads)/nReads^2
Xscg[,i] <- c(group.parameter[[i]][1:K], N_1Cj)
}

pi0 <- (Xscg[1:K,] %*% as.matrix(1/Xscg[K+1,]))/sum(1/Xscg[K+1,])
Xmcupo <- sum((1/Xscg[K+1,]) * apply((Xscg[1:K,]-pi0 %*% matrix(1, 1, numGrps))^2, 2, function(x){sum(x/pi0)}))

return(Xmcupo)
}
