library(e1071)

k_means=function(x,k){
  
  set.seed(42)
  x=as.matrix(x)
  res=kmeans(x,k)
  tss=res$tot.withinss
  
  for (i in 1:30){
    set.seed(42)
    res1=kmeans(x,k)
    tss1=res1$tot.withinss
    if (tss<tss1){res=res}else{res=res1}
    tss=res$tot.withinss
    
  }
  #return(res$centers)
  structure(list(centers=res$centers,cluster=res$cluster))
  
  
}


