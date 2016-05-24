iteration0 <-
function(N, p, m, d, h0, vc.method, poly.index, scale.h, lim.binning, nbins_, plmm, ...){
  arg<-list(...)
  ni<-plmm$ni
### E[y|T]
  Ey_t<-sm.regression(y=plmm$y, x=plmm$T_mat, h=h0[,1], eval.grid=F, eval.points=plmm$T_mat, display="none")$estimate
  
### E[X|T]
  EX_t<-NULL
  for(i in 1:p){
    EX.p_t<-sm.regression(y=plmm$X[,i], x=plmm$T_mat, h=h0[,i+1], eval.grid=F, eval.points=plmm$T_mat, display="none")$estimate
    EX_t<-cbind(EX_t, EX.p_t)
  }

### Tilde Matrix ###
  y_tilde<-plmm$y-Ey_t
  X_tilde<-plmm$X-EX_t
  plmm$y_tilde0<-split(y_tilde, plmm$cls)
  plmm$X_tilde0<-split(X_tilde, plmm$cls)
  
### VC estimation ###
  if(vc.method=="SA" || vc.method=="FC"){  
    Py<-ave(y_tilde, plmm$cls)
    PX<-matrix(NA, nrow=N, ncol=p)
    for(i in 1:p){
      PX[,i]<-ave(X_tilde[,i], plmm$cls)
    }  
    
### Estimate var_e 
    Qy<-y_tilde-Py# sorted using the order of cls
    QX<-X_tilde-PX# sorted using the order of cls
  
    invXQX<-solve(crossprod(X_tilde, QX)) # = solve(t(X_tilde)%*%QX)
    XQy<-crossprod(X_tilde, Qy)
    e_w<-c(y_tilde-X_tilde%*%invXQX%*%XQy) 
# invXQX%*%XQy = the same beta for cls-sorted and original => e_w is the same for cls-sorted and original   
    #e_w_bar=tapply(e_w, plmm$sort, mean)
    e_w_bar<-tapply(e_w, plmm$cls, mean)
    eQe<-e_w%*%e_w - sum(ni*e_w_bar*e_w_bar)#=ee-ePe
# e_w needs be unsorted. In an unbalanced case, the ordering of ni and e_w_bar must coincide.
    var_e<-c(eQe/(N-m-p))
     
    if(vc.method=="SA"){
### SA var_u
#      invXPX=solve(tX%*%PX)  
      #      XPy=tX%*%Py
      invXPX<-solve(crossprod(X_tilde, PX))
      XPy<-crossprod(X_tilde, Py)
      v_b<-y_tilde-X_tilde%*%invXPX%*%XPy
      #v_b_bar=tapply(v_b, plmm$sort, mean)
      v_b_bar<-tapply(v_b, plmm$cls, mean)
      vPv<-sum(ni*v_b_bar*v_b_bar)

      sum_xixi<-sapply(as.list(plmm$clsLev), SA_xixi, p=p, plmm=plmm)
      sum_xixi<-matrix(unlist(sum_xixi), nrow=p*p)
      sum_xixi<-matrix(rowSums(sum_xixi), ncol=p)

      var_u<-max( (vPv-(m-p)*var_e)/(N-sum(diag(invXPX%*%sum_xixi))), 0)

### FC var_u 
    }else if(vc.method=="FC"){
      v<-lm(y_tilde~X_tilde-1)$residuals
      vv<-v%*%v
      #XX<-t(X_tilde)%*%X_tilde # (p by p)
      XX<-crossprod(X_tilde) # (p by p)
      sum_nixx<-sapply(as.list(plmm$clsLev), FC_xixi, p=p, plmm=plmm)
      sum_nixx<-matrix(unlist(sum_nixx), nrow=p*p)
      sum_nixx<-matrix(rowSums(sum_nixx), ncol=p)
  
      if(p==1){
        var_u<-c( (vv-(N-p)*var_e)/(N-sum_nixx/XX) )
      }else{
        var_u<-c( (vv-(N-p)*var_e)/(N-sum(diag(solve(XX)%*%sum_nixx))) )  
      }
    var_u<-max(var_u, 0)
    }
  }
 
  if(vc.method=="FChetero"){      
    alpha<-plmm$alpha
    ay<-y_tilde/alpha
    aX<-X_tilde/alpha
    v_bar<-lm(ay~aX-1)$residuals
    vv<-v_bar%*%v_bar
      
    eQe_eta_Cfc<-FC_hetero(N=N, p=p, m=m, plmm=plmm, invarIND=NULL)
    eQe<-eQe_eta_Cfc$eQe
    eta_Cfc<-eQe_eta_Cfc$eta_Cfc

    var_e<-eQe/(N-m-p)
    var_u<-max( (vv-(N-p)*var_e)/(eta_Cfc), 0) 
    
  }else if(vc.method=="SAhetero"){
    eQe_vPv_eta_Csa<-SA_hetero(N=N, p=p, m=m, plmm=plmm, invarIND=NULL)
    eQe<-eQe_vPv_eta_Csa$eQe
    vPv<-eQe_vPv_eta_Csa$vPv
    eta_Csa<-eQe_vPv_eta_Csa$eta_Csa

    var_e<-eQe/(N-m-p)
    var_u<-max( (vPv-(m-p)*var_e)/eta_Csa, 0)
  }
#  VC0<-list(var_u=var_u, var_e=var_e)
  VC<-c(var_u=var_u, var_e=var_e)
  
### beta estimation ###
  beta<-est_beta(VC=VC, vc.method=vc.method, p=p, plmm=plmm)

### gamma Estimation ###
  yXb<-plmm$y-plmm$X%*%beta

  h<-h.select(y=yXb, x=plmm$T_mat, method="cv", poly.index=poly.index, ...)
  h<-h*scale.h 
  gamma<-sm.regression(y=yXb, x=plmm$T_mat, h=h, eval.grid=F, eval.points=plmm$T_mat, display="none", poly.index=poly.index)$estimate
  
  if(N<lim.binning){
    df_gamma<-cal_df_gamma(N=N, d=d, h=h, plmm=plmm, poly.index=poly.index)
  }else{ # binning
    if(d==1){
      df_gamma<-cal_df_gamma1(yXb=yXb, poly.index=poly.index, h=h, nbins=nbins_, plmm=plmm)
    }else if(d==2){
      df_gamma<-cal_df_gamma2(yXb=yXb, poly.index=poly.index, h=h, nbins=nbins_, plmm=plmm)
    }
  }
  
  return(list(beta=beta, VC=VC, gamma=gamma, h=h, df_gamma=df_gamma))
}
