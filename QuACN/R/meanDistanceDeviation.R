meanDistanceDeviation <- function(g, dist=NULL){
   # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  nV <- numNodes(g)
  
  dCM = rowSums(dist)
  wien = wiener(g,dist)
  
  Jii = sapply(1:nV, function(i){
    abs((dCM[i])-(2*wien/nV))
  })
    
  Ji <- abs(dCM-(2*wien/nV))
  return(1/nV * sum(unlist(Ji)))
}
