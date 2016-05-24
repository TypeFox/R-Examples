makedmat <-
function(nnod){
  ##creates designmatrix D with all possible assignments of the terminal nodes to the three  partition classes
  ##nnod = I = number of terminal nodes after a split 
  ##dmat=K * I matrix: K=number of possible assignments; 
  rmat<-3^(nnod)
  #rmat is total number of rows
  library(stats)
  dmat<-matrix(unlist(lapply(1:nnod,function(jj,rmat){as.double(gl(3,rmat/(3^jj),rmat))},rmat=rmat)),ncol=nnod,nrow=rmat)
  return(dmat)}
