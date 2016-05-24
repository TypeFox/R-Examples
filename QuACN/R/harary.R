harary <- function(g, dist=NULL){
  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  return ( sum(1/(dist[dist!=0]))/2)
}
