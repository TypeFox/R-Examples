LiStengos <-
function(v, var_u, LS.bws, poly.index, nonpar.bws, N, p, d, nbins_, LS.poly.index, trim, scale.h, lim.binning, plmm, ...){

### Li and Stengos Var function estimation ###
  #LS_h_var.ij=.var.ij(v=v, hetMat=hetMat, LS.bws=LS.bws, LS.poly.index=LS.poly.index, N=N, d=d,...)
  LS_h_var.ij<-.var.ij(v=v, LS.bws=LS.bws, LS.poly.index=LS.poly.index, N=N, d=d, plmm=plmm, ...)

  nu2_trim<-LS_h_var.ij$var.ij-var_u
  nu2_trim[nu2_trim < trim]<-trim # trim almost negative nu2 values
  #ifelse(nu2_trim < trim, trim, nu2_trim) # probably slower than the above
  plmm$nu2<-split(nu2_trim, plmm$cls)

### beta estimation ###
  XVXXVy<-sapply(as.list(plmm$clsLev), XVXXVy_LS, var_u=var_u, p, plmm=plmm)
  XVX<-XVXXVy[rownames(XVXXVy)=="XVX"]
  XVX<-matrix(unlist(XVX), nrow=p*p)
  XVX<-matrix(rowSums(XVX), ncol=p)

  XVy<-XVXXVy[rownames(XVXXVy)=="XVy"]
  XVy<-matrix(unlist(XVy), nrow=p)
  XVy<-rowSums(XVy)

  betaLS<-solve(XVX)%*%XVy

### gamma estimation ###
  yXb<-plmm$y-plmm$X%*%betaLS

  if(nonpar.bws=="h.select"){
    h<-h.select(y=yXb, x=plmm$T_mat, method="cv", poly.index=poly.index, ...)
  }else if(nonpar.bws=="hcv"){
    if(d==1){
      h<-hcv(y=yXb, x=as.vector(plmm$T_mat), poly.index=poly.index, ...)
    }else if(d==2){
      h<-hcv(y=yXb, x=plmm$T_mat, poly.index=poly.index, ...)
    } 
        
  }else if(nonpar.bws=="GCV"){
    sd_T<-apply(plmm$T_mat,2,sd)
    if(d==1){
      h_start<-apply(plmm$T_mat,2,sd)/2
    }else{ h_start<-0.5 }
      log_h<-optimize(f=GCV, interval=log(c(h_start/8, h_start*4)), yXb=yXb, d=d, N=N, poly.index=poly.index, plmm=plmm, sd_T=sd_T)
#  log_h=optimize(f=LS_GCV, interval=log(c(h_start/8, h_start*4)), yXb=yXb, var_ij=cal_nu2_$var.ij, var_u=var_u, d=d, N=N, poly.index=poly.index, nonpar.bws=nonpar.bws, plmm=plmm)#old LS_GCV using sm.regression
    h<-exp(log_h$minimum)
    if(d==2){h<-h*sd_T}      

  }else if(nonpar.bws=="GCV.c"){
    sd_T<-apply(plmm$T_mat,2,sd)
    if(d==1){
      h_start<-sd_T/2
    }else{ h_start<-0.5 }
    plmm$sd_ij0<-split(sqrt(nu2_trim+var_u), plmm$cls)
#browser()
    plmm$R.i_hetero<-lapply(plmm$clsLev, .R.i_hetero, var_u=var_u, plmm=plmm)# used in LScal_sRs.i and LScal_vech_R
    
    plmm$GCV_c<-NULL
    plmm$df_gamma<-NULL

    log_h<-optimize(f=.LS_GCV_c, interval=log(c(h_start/8, h_start*4)), yXb=yXb, var_u=var_u, d=d, N=N, poly.index=poly.index, plmm=plmm, sd_T=sd_T)
    #log_h=optimize(f=LS_GCV_c, interval=log(c(h_start/8, h_start*4)), yXb=yXb, var_ij=cal_nu2_$var.ij, var_u=var_u, d=d, N=N, poly.index=poly.index, nonpar.bws=nonpar.bws, plmm=plmm)#old LS_GCV_c using sm.regression

    h<-exp(log_h$minimum)
    
    #plmm$GCV_c=unlist(plmm$GCV_c)
    #plmm$df_gamma=unlist(plmm$df_gamma)
    df_gamma<-plmm$df_gamma[which(plmm$GCV_c==min(plmm$GCV_c))][1]
# df_gamma may have more than one same values if some GCV_c values are the same 
    if(d==2){h<-h*sd_T}
  }

  h<-h*scale.h
  gamma<-sm.regression(y=yXb, x=plmm$T_mat, h=h, eval.grid=F, eval.points=plmm$T_mat, display="none", poly.index=poly.index)$estimate

### check the estimates ###
#sm.regression(y=yXb, x=plmm$T_mat, h=h, poly.index=poly.index)
#lines(y=gamma_hat[order(plmm$T_mat)], x=sort(plmm$T_mat), col=4,lwd=2)

### Calculate df_gamma
#  if(alt.df==FALSE){
  if(nonpar.bws!="GCV.c"){
    if(N<lim.binning){ # without binning
      df_gamma<-cal_df_gamma(N=N, d=d, h=h, plmm=plmm, poly.index=poly.index)
    }else{ # binning
      if(d==1){
        df_gamma<-cal_df_gamma1(yXb=yXb, poly.index=poly.index, h=h, nbins=nbins_, plmm=plmm)
      }else if(d==2){
        df_gamma<-cal_df_gamma2(yXb=yXb, poly.index=poly.index, h=h, nbins=nbins_, plmm=plmm)
      }
    }
#  }else{# Carmack df_gamma
#    df_gamma=LScal_tr(var_ij=cal_nu2_$var.ij, var_u=var_u, d=d, h=h, N=N, nonpar.bws="GCV.c", plmm=plmm)
  }
   
  return(list(beta=betaLS, gamma=gamma, nu2=nu2_trim, h=h, LS_h=LS_h_var.ij$LS_h, df_gamma=df_gamma))
}
