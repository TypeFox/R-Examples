matrixrank <- function(x, tolerance=1e-8){
  #this function returns the rank of the matrix
  
  xsvd<-svd(x);
  xd<-xsvd$d;
  #rm(x.svd);
  xtrued<-xd[abs(xd)>tolerance]
  result<-length(xtrued)
  return(result)
}
