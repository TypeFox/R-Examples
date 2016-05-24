rGibbs <-
function(gstar, tau, type, numGraphs=1){
if(missing(gstar) || missing(type) || missing(tau))
stop("gstar, tau, prob, and/or type is missing.")

if(numGraphs <= 0)
stop("numGraphs must be an integer greater than 0.")

nodes <- getNumNodes(gstar, type)
edges <- getNumEdges(nodes, type)
gstar <- as.vector(as.matrix(gstar))

if(tolower(type) == "adjmatrix")
gstar <- full2lt(gstar)

#Calculate the prob for a graph k-distance away and then sample distances from that
possEdges <- 0:edges
distProb <- exp(lgamma(edges+1) - lgamma(edges-possEdges+1) - lgamma(possEdges+1) - tau*possEdges-edges * log(1+exp(-tau)))
gdists <- sample(possEdges, numGraphs, replace=TRUE, prob=distProb)

#Take our starting graph and randomly flip-flop k nodes to get a graph k-distance away
genData <- matrix(gstar, length(gstar), numGraphs)
for(i in 1:numGraphs){ 
nodesChange <- sample(1:length(gstar), gdists[i], replace=FALSE)
genData[nodesChange, i] <- 1*!genData[nodesChange, i]
}

if(tolower(type) == "adjmatrix")
genData <- apply(genData, 2, function(x){lt2full(x)})

return(as.data.frame(genData)) 
}
