
lglk.nb=function(n,r,v){
  lglk=lgamma(n+1/v)-lgamma(n+1)-lgamma(1/v)-log(1+r*v)/v-n*log(1+1/r/v)
  return(rowSums(lglk))
}
lglk.ps=function(n,r){
  lglk=-r+n*log(r)-lgamma(n+1)
  return(rowSums(lglk))
}



lglk.ps.1c=function(n,s,t,c){
  sc=sweep(s,2,c[t],"+")
  m=exp(sc)*rowSums(n+1e-10)/rowSums(exp(sc))
  lglk=-m+n*log(m)-lgamma(n+1)
  lglk=rowSums(lglk)
  return(lglk)
}
lglk.ps.c=function(n,s,t,C){
  if(is.vector(C)) C=matrix(C,nrow=1)
  nK=nrow(C)
  nG=nrow(n)
  lglk=matrix(0,nG,nK)
  for(k in 1:nK){
    lglk[,k]=lglk.ps.1c(n=n,s=s,t,c=C[k,])
  }
  return(lglk)
}




lglk.nb.1c=function(n,s,t,v,c){
  sc=sweep(s,2,c[t],"+")
  m=est.nb.mu.mle.one(n,exp(sc),v)
  m=exp(sc)*m
  lglk=rowSums(lgamma(n+1/v)-lgamma(n+1)-lgamma(1/v)-log(1+m*v)/v-log(1+1/m/v)*n)
  return(lglk)
}
lglk.nb.c=function(n,s,t,v,C){
  if(is.vector(C)) C=matrix(C,nrow=1)
  nK=nrow(C)
  nG=nrow(n)
  lglk=matrix(0,nG,nK)
  for(k in 1:nK){
    lglk[,k]=lglk.nb.1c(n=n,s=s,t,v,c=C[k,])
  }
  return(lglk)
}


