cal_sRs.i <-
function(clsLevel, plmm){
#cal_sRs.i=function(clsLevel, VC, plmm){
#  var_e=as.numeric(VC$var_e)
#  var_u=as.numeric(VC$var_u)
  
  i<-(plmm$clsLev==clsLevel)
#  ni=plmm$ni[i]
  s.k<-plmm$s.k0[i][[1]]
#  Ri=matrix(var_u/(var_u+var_e), ncol=ni, nrow=ni)
#  diag(Ri)=1
  
#  sRs.i=s.k%*%Ri%*%s.k
  sRs.i<-s.k%*%plmm$R.i[i][[1]]%*%s.k
  return(sRs.i)
}
