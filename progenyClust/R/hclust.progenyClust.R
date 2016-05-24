hclust.progenyClust<-
  function(x,k,h.method='ward.D2',dist='euclidean',p=2,...){
    d<-dist(x,method=dist,p)
    h<-hclust(d,h.method,...)
    cluster<-cutree(h,k)
    output=list(cluster=cluster,tree=h,dist=d)
    return(output)
  }