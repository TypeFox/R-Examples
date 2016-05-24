SA_xixi <-
function(clsLev, p, plmm){
  i<-(plmm$clsLev==clsLev)          
  X_tilde.i<-matrix(plmm$X_tilde0[i][[1]], ncol=p)  
  xi_tilde.i.<-colSums(X_tilde.i)
  xixi<-xi_tilde.i.%*%t(xi_tilde.i.)
  return(xixi)
}
