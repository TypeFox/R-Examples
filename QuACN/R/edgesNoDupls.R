.edgesNoDupls <- function(g){
  am <- adjacencyMatrix(g)
  am[lower.tri(am)]<-0
  edges <- apply(am,1,function(x){
    which(x!=0)
  })
  return(edges)
}
