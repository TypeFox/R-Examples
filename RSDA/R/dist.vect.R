dist.vect <-
function(vector1,vector2) {
  dist<-0
  if( ncol(vector1) == ncol(vector2) ){
    dist<-norm.vect(vector1-vector2)
  }
  return(dist) 
}
