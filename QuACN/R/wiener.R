wiener <- function(g, dist=NULL){
  #print("Wiener...")
  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  return(sum(dist)/2)
}
