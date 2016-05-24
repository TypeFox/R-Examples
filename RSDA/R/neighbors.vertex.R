neighbors.vertex <-
function(vertex,Matrix,num.neig) {
  dist<-dist.vect.matrix(t(vertex),Matrix)
  index<-order(dist,decreasing = FALSE)
  neighbors<-Matrix[index[1:num.neig],]
  return( list(neighbors = neighbors ,
               order = index ,
               distance = dist))
}
