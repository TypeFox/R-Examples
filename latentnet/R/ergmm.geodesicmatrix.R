ergmm.geodesicmatrix<-function(model){
  Ym<-model[["Ym"]]
  # For the purpose of geodesic distance, dichotomize the network about its mean.
  Ym<-Ym>mean(Ym,na.rm=TRUE)
  Ym[is.na(Ym)]<-0
  Ym <- Ym|t(Ym)
  geodist(Ym, count.paths=FALSE, inf.replace=nrow(Ym))$gdist
}
