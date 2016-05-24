transform_SR <-
function(clsLev, var_u, var_e, p, plmm){
  i<-(plmm$clsLev==clsLev)
  ni<-plmm$ni[i]
  aInv<-plmm$aInv0[i][[1]]
  eta<-sum(aInv^2)
  
  ay<-aInv*plmm$y_tilde0[i][[1]] 
  aX<-aInv*matrix(plmm$X_tilde0[i][[1]], ncol=p)
#  Q=diag(ni)-(1-sqrt(var_e/(var_e+eta*var_u)))/eta*aInv%*%t(aInv) 
  Q<-(1-sqrt(var_e/(var_e+eta*var_u)))/eta*aInv%*%t(aInv)   
  y<-ay-Q%*%ay
  X<-aX-Q%*%aX
  
  #XX=t(X)%*%X
  #Xy=t(X)%*%y
  XX<-crossprod(X, X)
  Xy<-crossprod(X, y)
  return(list(XX=XX, Xy=Xy))
}
