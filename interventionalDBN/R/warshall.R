warshall<-function(M) {
  d<-sqrt(length(M))
  path<-M
  for (k in 1:d) {
    for (i in 1:d) {
      for (j in 1:d) {
        path[i,j]<-max(path[i,j],min(path[i,k],path[k,j]))
      }
    }
  }
  return(path)
}