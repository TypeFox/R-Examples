cv.mogavs <-
function(mogavs,nvar,data,y_ind,K=10,R=1,order=FALSE){
  
  ind<-which(rowSums(mogavs$archiveSet)==nvar)
  mods<-data.frame(matrix(mogavs$archiveSet[ind,],ncol=ncol(data)-1))
  modslm<-data.frame("archIndex"=numeric(0),"formula"=character(0),"CVerror"=numeric(0),"CVse"=numeric(0),stringsAsFactors = FALSE)
  for(i in 1:nrow(mods)){
    lmx<-mogavsToLinear(mods[i,],y_ind,data)
    cvx<-cvTools::repCV(lmx,K=K,R=R)
    modslm[i,"archIndex"]<-ind[i]
    modslm[i,"formula"]<-deparse(lmx$terms,width.cutoff=500)
    modslm[i,"CVerror"]<-cvx$cv
    modslm[i,"CVse"]<-cvx$se
  }
  if(order)modslm<-modslm[order(modslm$CVerror),]
  return(modslm)
}
