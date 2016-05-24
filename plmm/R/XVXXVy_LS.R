XVXXVy_LS <-
function(clsLevel, var_u, p, plmm){
  i<-(plmm$clsLev==clsLevel)
  nu2Inv<-1/plmm$nu2[i][[1]]

  #ViInv=diag(nu2Inv)-(nu2Inv)%*%t(nu2Inv)/(1/var_u+sum(nu2Inv))
  ViInv<-diag(nu2Inv)-tcrossprod(nu2Inv, nu2Inv)/(1/var_u+sum(nu2Inv))
### check the inverse ###
#  Vi=var_u*matrix(1,ncol=plmm$ni[i],nrow=plmm$ni[i])+diag(plmm$nu2[i][[1]])
#round(Vi%*%ViInv,5)
#########################
  X_tilde.i<-matrix(plmm$X_tilde0[i][[1]], ncol=p)
  #XV=t(X_tilde.i)%*%ViInv
  XV<-crossprod(X_tilde.i, ViInv)
  XVX<-XV%*%X_tilde.i
  XVy<-XV%*%plmm$y_tilde0[i][[1]]
  return(list(XVX=XVX, XVy=XVy))
}
