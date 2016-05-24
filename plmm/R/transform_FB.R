transform_FB <-
function(clsLev, var_u, var_e, p, plmm){
  i<-(plmm$clsLev==clsLev)
  ni<-plmm$ni[i] 
  y<-plmm$y_tilde0[i][[1]] 
  X<-matrix(plmm$X_tilde0[i][[1]], ncol=p)
  Q<-(1-sqrt(var_e/(var_e+ni*var_u)))*matrix(1/ni, ncol=ni, nrow=ni)   
  y<-y-Q%*%y
  X<-X-Q%*%X
  
  #XX=t(X)%*%X
  #Xy=t(X)%*%y
  XX<-crossprod(X,X)
  Xy<-crossprod(X,y)
  
  return(list(XX=XX, Xy=Xy))
}
