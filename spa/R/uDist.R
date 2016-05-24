"uDist" <-
function(g,k){
   n<-dim(g)[1]
   m=g
   d=g
   for(i in 2:k){
      m0<-apply(d>0,2,as.numeric)
      m<-m0%*%g
      diag(m)<-0
      m<-apply(m>0,2,as.numeric)
      m<-i*apply((m-m0)>0,2,as.numeric)
      d=m+d
    }
    for(i in 1:n)d[d[,i]==0,i]=Inf
    diag(d)=0
    return(d)
}

