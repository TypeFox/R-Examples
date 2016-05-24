cal_XQXXQy_axx <-
function(clsLev, p, plmm, invarIND){
#browser()
  i<-(plmm$clsLev==clsLev)
  aInv<-plmm$aInv0[i][[1]]
  eta<-sum(aInv^2)
  ni<-plmm$ni[i]
  Q<-diag(ni)-aInv%*%t(aInv)/eta
  
  X_bar<-aInv*matrix(plmm$X_tilde0[i][[1]], ncol=p)
  y_bar<-aInv*plmm$y_tilde0[i][[1]]
  
#  XQ=t(X_underbar)%*%Q
  axx<-t(X_bar)%*%X_bar
  
  if(any(invarIND)){ X_bar<-X_bar[,!invarIND] }
#    XQX=XQ%*%X_underbar
#    XQy=XQ%*%y_underbar
  XQ<-t(X_bar)%*%Q
  XQX<-XQ%*%X_bar  
  XQy<-XQ%*%y_bar  

  return(list(XQX=XQX, XQy=XQy, axx=axx))
}
