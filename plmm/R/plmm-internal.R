.gamma_S <-
function(k, yXb, d, h, N, plmm, poly.index){
  #T_star=plmm$T_mat-matrix(1, nrow=N, ncol=d)%*%diag(plmm$T_mat[k,], ncol=d) 
  T_mat<-plmm$T_mat
  if(d==1){
    T_star<-T_mat-T_mat[k,]
    K<-c(exp(-0.5*((T_star/h)^2)))
# make it a vector to compute KT=K*cbind(1, T_star) in s.k.r
  }else{
    T_star<-T_mat-rep(T_mat[k,], each=N)
    #K=exp(-0.5*((T_star%*%diag(1/h, ncol=d))^2))
    K<-exp(-0.5*((T_star*rep(1/h, each=N))^2))
    K<-apply(K, 1, prod)
  } 
  
  s.k<-.s.k(T_star=T_star, K=K, d=d, poly.index=poly.index)
#  if(poly.index==0){ # NW
#    T_star=matrix(1, nrow=N, ncol=1)
#  }else{T_star=cbind(1, T_star)} #Local Linear
#  TK=t(K*T_star)
#  TKTinv=solve(TK%*%T_star)
#  s.k=(TKTinv%*%TK)[1,] # k-th row of S
 
  S.kk<-s.k[k]# S[k,k]
  gamma.k<-s.k%*%yXb # k-th gamma
  return(list(gamma.k=gamma.k, S.kk=S.kk))
}
.gamma_SR_SRS <-
function(k, yXb, d, h, N, plmm, poly.index){
#  T_star=plmm$T_mat-matrix(1, nrow=N, ncol=d)%*%diag(plmm$T_mat[k,], ncol=d)
#  K_multi=exp(-0.5*((T_star%*%diag(1/h, ncol=d))^2)) 
  T_mat<-plmm$T_mat
  if(d==1){
    T_star<-T_mat-T_mat[k,]
    K<-c(exp(-0.5*((T_star/h)^2)))
## make it a vector to compute KT=K*cbind(1, T_star) in s.k.r
  }else{
    T_star<-T_mat-rep(T_mat[k,], each=N)
    #K=exp(-0.5*((T_star%*%diag(1/h, ncol=d))^2))
    K<-exp(-0.5*((T_star*rep(1/h, each=N))^2))
    K<-apply(K, 1, prod)
  } 
  
  s.k<-.s.k(T_star=T_star, K=K, d=d, poly.index=poly.index)
#  K=apply(K_multi, 1, prod)
#  if(poly.index==0){ # NW
#    T_star=matrix(1, nrow=N, ncol=1)
#  }else{T_star=cbind(1, T_star)} #Local Linear
#  TK=t(K*T_star)
#  TKTinv=solve(TK%*%T_star)
#  s.k=(TKTinv%*%TK)[1,] # k-th row of S

  gamma.k<-s.k%*%yXb # k-th gamma
  plmm$s.k0<-split(s.k, plmm$cls)
  
### tr(SRS)
#  sRs.k=sapply(as.list(plmm$clsLev), cal_sRs.i, VC=VC, plmm=plmm)#calculatin for sum(s.ki%*%R.i%*%s.ki)
  sRs.k<-sapply(as.list(plmm$clsLev), cal_sRs.i, plmm=plmm)
  sRs.k<-sum(sRs.k)

### For tr(SR)
  k.cls<-(plmm$clsLev==plmm$cls[k])# cls of obs k 
## trim down s.k to its relevant part 
  s.ki<-s.k[rep(k.cls, times=plmm$ni)]

  return(list(gamma.k=gamma.k, sRs.k=sRs.k, s.ki=s.ki))
}
.GCV_c <-
function(h, yXb, d, N, poly.index, plmm, sd_T){
  h<-exp(h)
  if(d==2){h<-h*sd_T}
#  gamma=sm.regression(y=yXb, x=plmm$T_mat, h=h, eval.grid=F, eval.points=plmm$T_mat, poly.index=poly.index, display="none")$estimate

  gamma_SR_SRS<-sapply(as.list(1:N), .gamma_SR_SRS, yXb=yXb, d=d, h=h, N=N, poly.index=poly.index, plmm=plmm)
  gamma<-unlist(gamma_SR_SRS[rownames(gamma_SR_SRS)=="gamma.k"])
  #gamma=unlist(gamma_SR_SRS[1,])
  trSRS<-sum(unlist(gamma_SR_SRS[rownames(gamma_SR_SRS)=="sRs.k"]))
  #trSRS=sum(unlist(gamma_SR_SRS[2,]))
  vech_s<-unlist(gamma_SR_SRS[rownames(gamma_SR_SRS)=="s.ki"])
  
  #vech_s=unlist(gamma_SR_SRS[3,])
#  vech_R=sapply(as.list(plmm$clsLev), cal_vech_R, var_u, var_e, plmm=plmm)
#  if(is.list(vech_R)){vech_R=unlist(vech_R)
#  }else{vech_R=c(vech_R)}
  vech_R<-unlist(plmm$R.i)
  trSR<-vech_s%*%vech_R
  #trSR=vech_s%*%plmm$vech_R

  DF<-2*trSR-trSRS
  GCV_c<-sum((yXb-gamma)*(yXb-gamma))/(N*(1-DF/N)^2)
  
#  plmm$df_gamma=list(unlist(plmm$df_gamma), DF)
#  plmm$GCV_c=list(unlist(plmm$GCV_c), GCV_c)
  plmm$df_gamma<-c(plmm$df_gamma, DF)
  plmm$GCV_c<-c(plmm$GCV_c, GCV_c)
#  plmm$h=list(unlist(plmm$h),log(h))
  
  return(c(GCV_c))  
#  DF=cal_tr(VC=VC, d=d, h=h, N=N, nonpar.bws=nonpar.bws, plmm=plmm, poly.index=poly.index)
}
.LS_GCV_c <-
function(h, yXb, var_u, d, N, poly.index, plmm, sd_T){
  h<-exp(h)
  if(d==2){h<-h*sd_T}  
#  gamma=sm.regression(y=yXb, x=plmm$T_mat, h=h, eval.grid=F, eval.points=plmm$T_mat, poly.index=poly.index, display="none")$estimate

  LSgamma_SR_SRS<-sapply(as.list(1:N), .LSgamma_SR_SRS, yXb, var_u=var_u, d=d, h=h, N=N, poly.index=poly.index, plmm=plmm)

  gamma<-unlist(LSgamma_SR_SRS[rownames(LSgamma_SR_SRS)=="gamma.k"])
  trSRS<-sum(unlist(LSgamma_SR_SRS[rownames(LSgamma_SR_SRS)=="sRs.k"]))
  vech_s<-unlist(LSgamma_SR_SRS[rownames(LSgamma_SR_SRS)=="s.ki"])
  
  vech_R<-unlist(plmm$R.i_hetero)
  trSR<-vech_s%*%vech_R
  #trSR=vech_s%*%plmm$vech_R
  DF<-2*trSR-trSRS   
  
  GCV_c<-sum((yXb-gamma)*(yXb-gamma))/(N*(1-DF/N)^2)
  
#  plmm$df_gamma<-list(unlist(plmm$df_gamma), DF)
#  plmm$GCV_c=list(unlist(plmm$GCV_c), GCV_c)
  plmm$df_gamma<-c(plmm$df_gamma, DF)
  plmm$GCV_c<-c(plmm$GCV_c, GCV_c)
  
  return(c(GCV_c))
}
.LSgamma_SR_SRS <-
function(k, yXb, var_u, d, h, N, poly.index, plmm){
  #T_star=plmm$T_mat-matrix(1, nrow=N, ncol=d)%*%diag(plmm$T_mat[k,], ncol=d)
  if(d==1){
    T_star<-plmm$T_mat-plmm$T_mat[k,]
    K<-c(exp(-0.5*((T_star/h)^2)))
## make it a vector to compute KT=K*cbind(1, T_star) in s.k.r
  }else{
    T_star<-plmm$T_mat-rep(plmm$T_mat[k,], each=N)
    #K=exp(-0.5*((T_star%*%diag(1/h, ncol=d))^2))
    K<-exp(-0.5*((T_star*rep(1/h, each=N))^2))
    K<-apply(K, 1, prod)
  } 
  
  s.k<-.s.k(T_star=T_star, K=K, d=d, poly.index=poly.index)
#  K_multi=exp(-0.5*((T_star%*%diag(1/h, ncol=d))^2)) 
#  K=apply(K_multi, 1, prod)
#  if(poly.index==0){ # NW
#    T_star=matrix(1, nrow=N, ncol=1)
#  }else{T_star=cbind(1, T_star)} #Local Linear
#  TK=t(K*T_star)
#  TKTinv=solve(TK%*%T_star)
#  s.k=(TKTinv%*%TK)[1,] # k-th row of S
  gamma.k<-s.k%*%yXb # k-th gamma

  plmm$s.k0<-split(s.k, plmm$cls)

  #  plmm$var_ij0=split(var_ij, plmm$cls)
  
### tr(SRS)
#  sRs.k=sapply(as.list(plmm$clsLev), LScal_sRs.i, var_u=var_u, plmm=plmm)#calculatin for sum(s.ki%*%R.i%*%s.ki)
  sRs.k<-sapply(as.list(plmm$clsLev), LScal_sRs.i, plmm=plmm)#calculatin for sum(s.ki%*%R.i%*%s.ki)
  sRs.k<-sum(sRs.k)

### For tr(SR)
  k.cls<-(plmm$clsLev==plmm$cls[k])# cls of obs k 
## trim down s.k to its relevant part 
  s.ki<-s.k[rep(k.cls,times=plmm$ni)]      

  return(list(gamma.k=gamma.k, sRs.k=sRs.k, s.ki=s.ki))
}
.R.i <-
function(clsLev, VC, plmm){
  i<-(plmm$clsLev==clsLev)
  ni<-plmm$ni[i]
  Ri<-matrix(VC[1]/(VC[1]+VC[2]), ncol=ni, nrow=ni)
  diag(Ri)<-1
  return(Ri)  
}
.R.i_hetero <-
function(clsLev, var_u, plmm){
  i<-(plmm$clsLev==clsLev)
  ni<-plmm$ni[i]
  Ri<-matrix(var_u, ncol=ni, nrow=ni)
  A<-diag(1/plmm$sd_ij0[i][[1]])
  Ri<-A%*%Ri%*%A  
  diag(Ri)<-1
  return(Ri)  
}
.s.k <-
function(T_star, K, d, poly.index){
  if(poly.index==1){# Local Linear
    if(d==1){
      TKT11<-sum(K)
      TKT12<-sum(K*T_star)  
      TKT22<-sum(K*T_star*T_star)
      TKTinv1<-c(TKT22, -TKT12)/(TKT11*TKT22-TKT12^2)
         
#      TK=rbind(1, T_star)*rep(K, each=2)
      KT<-K*cbind(1, T_star)
    }else if(d==2){      
      TKT11<-sum(K)
      TKT12<-sum(K*T_star[,1])
      TKT13<-sum(K*T_star[,2])
      TKT22<-sum(K*T_star[,1]*T_star[,1])
      TKT23<-sum(K*T_star[,1]*T_star[,2])
      TKT33<-sum(K*T_star[,2]*T_star[,2])
      det<-TKT11*TKT22*TKT33+2*TKT12*TKT23*TKT13-(TKT22*TKT13^2+TKT11*TKT23^2+TKT33*TKT12^2)
    
      TKTadj11<-TKT22*TKT33-TKT23^2
      TKTadj12<--(TKT12*TKT33-TKT13*TKT23)
      TKTadj13<-TKT12*TKT23-TKT13*TKT22
      TKTinv1<-c(TKTadj11, TKTadj12, TKTadj13)/det
     
#      TK=rbind(1, T_star[,1], T_star[,2])*rep(K, each=3) 
      KT<-K*cbind(1, T_star) 
    }
  }else{ # N.W.
    TKTinv1<-1/sum(K)
#    TK=K
    KT<-K # (N by 1)
  }
#  return(TKTinv1%*%TK) # = TKTinv%*%TK[1,] = S[k,] = s.k
  return(c(KT%*%TKTinv1))
## return a vector to compute gamma.k=s.k%*%yXb in GCV.r
}
.var.ij <-
function(v, LS.bws, LS.poly.index, N, d, plmm, ...){
  arg<-list(...)
  hetMat<-plmm$hetMat
  if(!is.null(arg$nbins)){nbins<-arg$nbins
  }else{nbins<-round(8*log(N)/d)}
   
### nu.ij ###
  v2<-v^2# for variance function estimation
  if(LS.bws=="hcv"){
    if(is.vector(hetMat)){
      LS_h<-hcv(y=v2, x=as.vector(hetMat), poly.index=LS.poly.index)    
    }else{
      LS_h<-hcv(y=v2, x=hetMat, poly.index=LS.poly.index)
    }    
  }else if(LS.bws=="h.select"){
    LS_h<-h.select(y=v2, x=hetMat, nbins=nbins, method="cv", poly.index=LS.poly.index, ...)
  }else if(LS.bws=="ROT"){
    if(is.vector(hetMat)){
      LS_h<-sd(hetMat)*N^(-1/(4+1))
    }else{
      dimH<-ncol(hetMat)
      LS_h<-apply(hetMat, 2, sd)*N^(-1/(4+dimH))
    }
  }
# <Experiment> Check the var.fun.bws performance
#sm.regression(y=v2, x=hetMat, h=h)

  var.ij<-sm.regression(y=v2, x=hetMat, h=LS_h, eval.grid=F, eval.points=hetMat, poly.index=LS.poly.index, display="none")$estimate
#sm.regression(y=v2, x=hetMat, h=h_LS, poly.index=0)

#  nu2=var.ij-var_u
#  return(list(nu2=nu2, LS_h=LS_h, var.ij=var.ij))
  return(list(LS_h=LS_h, var.ij=var.ij))
}
