cluster.size <-
function(id){
  clid<- unique(id)
  m<- length(unique(id))
  n<- rep(0,m)
  autotime<- rep(0,0)
  for(i in 1:m){
    n[i]<- length(which(id==clid[i]))
    autotime<- c(autotime,1:n[i])
  }
  id<- rep(1:m,n)
  return(list(m=m,n=n,id=id,autotime=autotime))
}
