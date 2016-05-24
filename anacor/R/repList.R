repList<-function(x,n) {
  z<-list()
  for (i in 1:n)
    z<-c(z,list(x))
  return(z)
}