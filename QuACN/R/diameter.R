diameter <- function(g, dist=NULL){
  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' has to be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))

  
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  return(max(dist))
}
