productOfRowSums <- function(g, dist=NULL, log=FALSE){
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  if(log==FALSE){
    return(prod(apply(dist,1,sum)))
  }else{
    return(log2(prod(apply(dist,1,sum))))
  }
}
