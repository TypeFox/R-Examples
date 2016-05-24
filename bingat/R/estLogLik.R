estLogLik <-
function(data, type, g, tau){
if(missing(data) || missing(type) || missing(g) || missing(tau))
stop("data, type, g, and/or tau is missing.")

nodes <- getNumNodes(data, type)
edges <- getNumEdges(nodes, type)
normConst <- 1+exp(-tau)    
###normConst <- ifelse(normConst==0, .Machine$double.xmin, normConst) #Adjust normConst if it is 0

if(!is.null(dim(g))){ #Check gstar is a single vector or matrix
distToGStar <- NULL
for(i in 1:ncol(g)){
distToGStar[i] <- calcDistance(data[,i], g[,i], type)}
}else{
distToGStar <- apply(data, 2, function(x, g, type) {calcDistance(x, g, type)}, g=g, type=type)
}

LogLik <- -edges * ncol(data) * log(normConst) - tau * sum(distToGStar)

return(LogLik)
}
