adjacencyMatrix <- function(g){
  # check if g is a graphNEL object
   if(class(g)[1]!="graphNEL"){
    stop("'g' has to be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  am <- as(g, "matrix")
  if(any(am>1)){
    warning("QuACN currently only supports unweighted graphs. The matrix will be unweighted")
    am <- ifelse(am>1,1,am)
  }
  return(am)
}
