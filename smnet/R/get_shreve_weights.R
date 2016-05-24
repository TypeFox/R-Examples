
get_shreve_weights<-function(adjacency, shreve.order){
  # convert shreve to weights
  n.segments<-nrow(adjacency)
  shreve.weights<-vector("numeric", length = nrow(adjacency))
  ij<-triplet(adjacency)$indices
  
  for(j in 1:n.segments){
    ij.ind<-which(ij[,2] == j)
    shreve.down<-shreve.order[j]
    if(length(ij.ind)>0){
      for(i in 1:length(ij.ind)){
        shreve.up<-shreve.order[ij[ij.ind[i],1]]
        shreve.weights[ij[ij.ind[i]]]<-shreve.up/shreve.down
      }
    }
  }
  shreve.weights[which(shreve.weights == 0)]<-1# this shouldn't happen and so there is an error....
  shreve.weights
}