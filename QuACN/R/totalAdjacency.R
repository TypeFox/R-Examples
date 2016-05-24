totalAdjacency <- function(g, am=NULL){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(am)){
    am <- adjacencyMatrix(g)
  }
  sum(am)/2
}
