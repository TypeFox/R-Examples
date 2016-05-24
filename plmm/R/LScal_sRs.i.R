LScal_sRs.i <-
function(clsLevel, plmm){  
#LScal_sRs.i=function(clsLevel, var_u, plmm){  
  i<-(plmm$clsLev==clsLevel)
#  ni=plmm$ni[i]
  s.k<-plmm$s.k0[i][[1]]
#  Ri=matrix(var_u, ncol=ni, nrow=ni)
#  A=diag(1/plmm$var_ij0[i][[1]])
#  A=diag(1/plmm$sd_ij0[i][[1]])
#  Ri=A%*%Ri%*%A  
#  diag(Ri)=1  
#  sRs.i=s.k%*%Ri%*%s.k
  sRs.i<-s.k%*%plmm$R.i_hetero[i][[1]]%*%s.k
  return(sRs.i)
}
