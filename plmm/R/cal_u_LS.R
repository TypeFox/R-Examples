cal_u_LS <-
function(clsLevel, var_u, LS, plmm){
  i<-(plmm$clsLev==clsLevel)
  ni<-plmm$ni[i]
  if(LS){ # wplmm
    nu2Inv<-plmm$nu2Inv0[i][[1]]
    Vinv<-diag(nu2Inv)-outer(nu2Inv, nu2Inv, "*")*var_u/(1+var_u*sum(nu2Inv))
  }else{ Vinv<-plmm$ViInv[i][[1]] } # plmm with hetero.prop
  #u_hat=var_u*rep(1, ni)%*%Vinv%*%(plmm$y_Xb_gamma0[i][[1]])
  u_hat<-var_u*rep(1, ni)%*%Vinv%*%(plmm$v0[i][[1]])
  return(u_hat)
}
