ClassifPaired<-function(Data,Tcla)
{
  ncrit<-length(Data@Paircomp)
  nsujet<-length(Data@Cons)
  dc<-as.dist(matrix(0,nsujet,nsujet))  
  for(k in 1:ncrit)
  {
    vk<-NULL
    for (h in 1:nsujet)
    {
      vk<-cbind(vk,simplify2array(lapply(Data@Paircomp[[k]][h],rowSums)))
    }
    vk<-t(vk)
    vk<-vk/rowSums(vk)
    dck<-dist(vk,method="euclidean",diag=FALSE,upper=FALSE)
    dc<-dc+dck
  }
  dc<-as.dist(dc)
  classif<-hclust(dc,method="ward.D2")
  appart<-cutree(classif,Tcla)
  return(appart)
  
}
