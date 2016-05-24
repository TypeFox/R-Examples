normalizedEdgeComplexity <- function(g,ita=NULL){
   if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(ita)){
    ita <-  totalAdjacency(g)
    }
  return(ita/(numNodes(g)^2))
}
