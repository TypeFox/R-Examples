make_weights<-function(adjacency, type = "unif"){
  adjacency<-as.matrix(adjacency)
  nrowAdj<-nrow(adjacency)
  weights<-vector("numeric", length = nrowAdj)
  if(type == "unif"){
    for(i in 1:nrowAdj){
      inds<-which(adjacency[,i] == 1)
      if(length(inds) == 1) weights[inds] <- 1
      if(length(inds) == 2){
        weights[inds[1]]<-runif(1)
        weights[inds[2]]<-1-weights[inds[1]]
      }
    }
  }
  return(weights)
}