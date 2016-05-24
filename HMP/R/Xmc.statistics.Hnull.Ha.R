Xmc.statistics.Hnull.Ha <-
function(K, pi0, group.parameter){
numGroups <- length(group.parameter)
group.parameter.estimated <- vector("list", numGroups)

for(x in 1:numGroups){
pi <- group.parameter[[x]][1:K]
theta <- group.parameter[[x]][K+1]
P <- length(group.parameter[[x]])
nreads.data <- as.matrix(group.parameter[[x]][(K+2):P])
data <- Dirichlet.multinomial(nreads.data, shape=pi*(1-theta)/theta)

pi.MoM <- colSums(data)/sum(colSums(data))
theta.MoM <- weirMoM(data, pi.MoM)

group.parameter.estimated[[x]] <- c(pi.MoM, theta.MoM, t(nreads.data))
}

Xmc <- Xmc.statistics(group.parameter.estimated, pi0)

return(Xmc)
}
