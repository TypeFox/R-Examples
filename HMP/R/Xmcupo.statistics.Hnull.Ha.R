Xmcupo.statistics.Hnull.Ha <-
function(K, group.parameter){
numGroups <- length(group.parameter)
group.parameter.estimated <- vector("list", numGroups)

for(x in 1:numGroups){
pi <- group.parameter[[x]][1:K]
theta <- group.parameter[[x]][K+1]
P <- length(group.parameter[[x]])
nreads.data <- as.matrix(group.parameter[[x]][(K+2):P])
data <- Dirichlet.multinomial(nreads.data, shape=pi*(1-theta)/theta)

totalsample <- sum(colSums(data))
pi.MoM <- colSums(data)/totalsample

q <- pi.MoM
r <- length(as.vector(which(q==0)))
rr <- length(q)-r
q[which(q!=0)] <- q[which(q!=0)]-r/(rr*2*(totalsample+1))
q[which(q==0)] <- 1/(2*(totalsample+1))
pi.MoMb <- q

theta.MoM <- weirMoM(data, pi.MoM)
group.parameter.estimated[[x]] <- c(pi.MoMb, theta.MoM, t(nreads.data))
}

Xmcupo <- Xmcupo.statistics(group.parameter.estimated, K)

return(Xmcupo)
}
