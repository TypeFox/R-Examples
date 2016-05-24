inv_Vi_hetero <-
function(clsLev, var_u, var_e, alpha, plmm){
  i<-(plmm$clsLev==clsLev)
  ni<-plmm$ni[i]
  aInv2<-plmm$aInv0[i][[1]]^2
  eta<-sum(aInv2)
  ViInv<-(diag(aInv2)-(var_u/(eta*var_u+var_e))*outer(aInv2, aInv2, "*"))/var_e
  return(ViInv)
}
