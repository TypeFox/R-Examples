cal_eQe_eta_Cfc <-
function(clsLev, beta_w, axxInv, p, plmm, invarIND){
#browser()
  i<-(plmm$clsLev==clsLev)
  aInv<-plmm$aInv0[i][[1]]
  eta<-sum(aInv^2)
  ni<-plmm$ni[i]
  Q<-diag(ni)-aInv%*%t(aInv)/eta

  X_bar<-aInv*matrix(plmm$X_tilde0[i][[1]], ncol=p)
  y_bar<-aInv*plmm$y_tilde0[i][[1]]

### Cfc
  Cfc<-aInv%*%(X_bar)%*%axxInv%*%t(X_bar)%*%aInv
  
### eQe
  if(any(invarIND)){ X_bar<-X_bar[,!invarIND] }
  e_w<-c(y_bar - X_bar%*%beta_w)
  # make it a vector, not a matrix so that next line work
#  SSE3=c(ei%*%solve(Qi[-ni,-ni])%*%ei)
  eQe<-c(e_w%*%Q%*%e_w)

  return(list(eQe=eQe, eta_Cfc=eta-Cfc))
}
