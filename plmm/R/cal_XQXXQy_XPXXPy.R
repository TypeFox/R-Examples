cal_XQXXQy_XPXXPy <-
function(clsLev, p, plmm, invarIND){
  i<-(plmm$clsLev==clsLev)
  aInv<-plmm$aInv0[i][[1]]
  eta<-sum(aInv^2)
  ni<-plmm$ni[i]

  X_bar<-aInv*matrix(plmm$X_tilde0[i][[1]], ncol=p)
  y_bar<-aInv*plmm$y_tilde0[i][[1]]
  
  P<-aInv%*%t(aInv)/eta
  Q<-diag(ni)-P
  XP<-t(X_bar)%*%P
  XPX<-XP%*%X_bar
  XPy<-XP%*%y_bar

  if(any(invarIND)){ X_bar<-X_bar[,!invarIND] }
  XQ<-t(X_bar)%*%Q
  XQX<-XQ%*%X_bar
  XQy<-XQ%*%y_bar
  
  return(list(XQX=XQX, XQy=XQy, XPX=XPX, XPy=XPy))
}
