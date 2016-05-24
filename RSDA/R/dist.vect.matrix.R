dist.vect.matrix <-
function(vector,Matrix) {
  num.col.vect<-ncol(vector)
  num.col.matrix<-ncol(Matrix)
  num.row.matrix<-nrow(Matrix)
  dist.matrix<-rep(-1,num.row.matrix)
  if(num.col.vect == num.col.matrix){
    for(i in 1:num.row.matrix){
      w<-as.matrix(Matrix[i,])
      dist.matrix[i]<-dist.vect(vector,t(w))
    }
  }
  return(dist.matrix)
}
