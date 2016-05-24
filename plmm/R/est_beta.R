est_beta <-
function(VC, vc.method, p, plmm){
#  var_e=as.numeric(VC$var_e)
#  var_u=as.numeric(VC$var_u)
  
  if(vc.method=="FC" || vc.method=="SA"){
    #XXXy=sapply(as.list(plmm$clsLev), transform_FB, var_u=var_u, var_e=var_e, p=p, plmm=plmm)
    XXXy<-sapply(as.list(plmm$clsLev), transform_FB, var_u=VC[1], var_e=VC[2], p=p, plmm=plmm)
  }else if(vc.method=="FChetero"|| vc.method=="SAhetero"){
#    XXXy=sapply(as.list(plmm$clsLev), transform_SR, var_u=var_u, var_e=var_e, p=p, plmm=plmm)
    XXXy<-sapply(as.list(plmm$clsLev), transform_SR, var_u=VC[1], var_e=VC[2], p=p, plmm=plmm)
  }
#  XXXy=list of 2m elements (m elements with "XX" and m elements with "Xy")
  
  XX<-XXXy[rownames(XXXy)=="XX"]
  XX<-matrix(unlist(XX), nrow=p*p)
  XX<-matrix(rowSums(XX), ncol=p)

  Xy<-XXXy[rownames(XXXy)=="Xy"]
  Xy<-matrix(unlist(Xy), nrow=p)
  Xy<-rowSums(Xy)
  
  beta<-solve(XX)%*%Xy
#  colnames(beta)="Estimate"
#  rownames(beta)=plmm$xName
  return(beta)
}
