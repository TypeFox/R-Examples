floyd<-function(x,verbose=FALSE){
  if(!is.matrix(x) & !is.data.frame(x))
    stop("x must be a matrix or a data.frame")
  m<-dim(x)[1]
  n<-m*m
  d=as.vector(as.matrix(x))
  baser=1:m
  basec=m*(0:(m-1))
  t1<-proc.time()
  ind1<-sort(rep(1:m,m))
  ind2<-rep(1:m,m)
  for(k in 1:m){
    m2<-d[baser+m*(k-1)][ind1]+d[baser+m*(k-1)][ind2]
    ind<-m2<d
    d[ind]<-m2[ind]
    if(verbose)
      cat("k=",k,".... time=",(proc.time()-t1)/60, "\n")
  } 
  return(matrix(d,m,m))
}
