getLoglikeMixture <-
function(data, mixture, numConst){ 
if(missing(data) || missing(mixture))
stop("data and/or mixture is missing.")

type <- mixture$type
numGroups <- mixture$numgroups
numGraphs <- ncol(data)
nodes <- getNumNodes(data, type)
edges <- getNumEdges(nodes, type)

#Calculate distances to the gstars
distances <- matrix(NA, ncol(data), numGroups)
for(j in 1:numGroups)
distances[,j] <- apply(data, 2, function(x){calcDistance(x, mixture$gstars[[j]], type)})

t <- matrix(NA, ncol(data), numGroups)
for(i in 1:numGraphs)
t[i,] <- log(mixture$weights) - edges*log((1+exp(-mixture$taus))) - mixture$taus*distances[i,]

LL <- sum(apply(t, 1, matrixStats::logSumExp))

if(missing(numConst))
numConst <- edges/nodes * 10 * log(numGraphs) #(2+edges) * log(numGraphs)
BIC <- -2*LL + numGroups * numConst

return(list(ll=LL, bic=BIC))
}
