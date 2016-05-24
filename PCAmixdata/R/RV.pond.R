RV.pond<-function(liste.mat,ponde){
  Lg<-Lg.pond(liste.mat,ponde)
  RV<-matrix(NA,ncol=ncol(Lg),nrow=nrow(Lg))
  for (i in 1: nrow(Lg)){
    for (j in 1:ncol(Lg)){
      RV[i,j]<-Lg[i,j]/sqrt(Lg[i,i]*Lg[j,j])
    }
  }
  rownames(RV)<-colnames(RV)<-rownames(Lg)
  return(RV)
}