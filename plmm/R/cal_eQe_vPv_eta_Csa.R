cal_eQe_vPv_eta_Csa <-
function(clsLev, beta_w, beta_b, p, invXPX, plmm, invarIND){
#browser()
  i<-(plmm$clsLev==clsLev)
  aInv<-plmm$aInv0[i][[1]]
  eta<-sum(aInv^2)
  ni<-plmm$ni[i]

  P<-aInv%*%t(aInv)/eta
  Q<-diag(ni)-P
  X_bar<-aInv*matrix(plmm$X_tilde0[i][[1]], ncol=p)
  y_bar<-aInv*plmm$y_tilde0[i][[1]]

### vPv
  v<-c(y_bar-X_bar%*%beta_b)
  vPv<-c(v%*%P%*%v)

### Csa
  Csa<-aInv%*%(X_bar)%*%invXPX%*%t(X_bar)%*%aInv  

### eQe
  if(any(invarIND)){ X_bar<-X_bar[,!invarIND] }
  e_w<-c(y_bar-X_bar%*%beta_w)
  eQe<-c(e_w%*%Q%*%e_w)

  return(list(eQe=eQe, vPv=vPv, eta_Csa=eta-Csa))
}
