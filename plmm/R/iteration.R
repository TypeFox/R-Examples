iteration <-
function(beta, VC, gamma, h, df_gamma, vc.method, iter, nonpar.bws, poly.index, N, m, p, d, epsilon, scale.h, lim.binning, nbins_, plmm, BS, ...){
  arg<-list(...)
  VC_iter<-matrix(NA, ncol=2, nrow=iter+1)
  colnames(VC_iter)<-c("var_u", "var_e")
  VC_iter[1,]<-VC
  beta_iter<-matrix(NA, ncol=p, nrow=iter+1)
  beta_iter[1,]<-beta  
  ni<-plmm$ni
  cls<-plmm$cls
  T_mat<-plmm$T_mat
  X_tilde<-plmm$X # Re-define X_tilde
  plmm$X_tilde0<-split(X_tilde, cls) 
# Re-define plmm$X_tilde0 for FC_hetero, SA_hetero, FC_xixi, SA_xixi
  
#### Preparation for iteration() ####
  if(!BS){
    PX<-matrix(NA, nrow=N, ncol=p)
    for(i in 1:p){
      PX[,i]<-ave(X_tilde[,i], cls)
    }
  
    QX<-X_tilde-PX # plmm$X = X_tilde in the iteration
  #invarIND=apply(QX, 2, function(x){identical(all.equal(0, sum(x%*%x)), T)})# invarIND can be used for hetero case
    invarIND<-apply(QX, 2, function(x){identical(all.equal(0, sum(x*x)), T)})
    if(any(invarIND)){ QX<-QX[,!invarIND] }
    ncolQX<-ncol(QX)
  
    if(vc.method=="SA" || vc.method=="FC"){
      if(any(invarIND)){# invariant variables
        invXQX<-solve(crossprod(X_tilde[,!invarIND], QX))
      }else{ invXQX<-solve(crossprod(X_tilde, QX)) }
    
      if(vc.method=="SA"){
        invXPX<-solve(crossprod(X_tilde, PX))
        
        sum_xixi<-sapply(as.list(plmm$clsLev), SA_xixi, p=p, plmm=plmm)
        sum_xixi<-matrix(unlist(sum_xixi), nrow=p*p)
        sum_xixi<-matrix(rowSums(sum_xixi), ncol=p)
      
      }else if(vc.method=="FC"){      
        XX<-crossprod(X_tilde) # (p by p)
      
        sum_nixx<-sapply(as.list(plmm$clsLev), FC_xixi, p=p, plmm=plmm)
        sum_nixx<-matrix(unlist(sum_nixx), nrow=p*p)
        sum_nixx<-matrix(rowSums(sum_nixx), ncol=p)
      }
    }
  }else{
    invarIND<-plmm$invarIND
    invXQX<-plmm$invXQX
    invXPX<-plmm$invXPX
    ncolQX<-plmm$ncolQX
    XX<-plmm$XX
    sum_xixi<-plmm$sum_xixi
    sum_nixx<-plmm$sum_nixx
  }
    
################ Iteration ################
  r<-0
  while(iter > r){ 
    y_tilde<-(plmm$y-gamma) # Re-define and update y_tilde
    plmm$y_tilde0<-split(y_tilde, cls)
    
### VC estimation ###
    if(vc.method=="SA" || vc.method=="FC"){
      Py<-ave(y_tilde, cls)
  
### Estimate var_e 
      Qy<-y_tilde-Py
      if(any(invarIND)){# there are invariant variables
        XQy<-crossprod(X_tilde[,!invarIND], Qy)
        e_w<-c(y_tilde-X_tilde[,!invarIND]%*%invXQX%*%XQy)   
      }else{
        XQy<-crossprod(X_tilde, Qy)
        e_w<-c(y_tilde-X_tilde%*%invXQX%*%XQy)   
      }

      e_w_bar<-tapply(e_w, cls, mean)
      eQe<-e_w%*%e_w - sum(ni*e_w_bar*e_w_bar)#=ee-ePe
      var_e<-c(eQe/(N-m-ncolQX-df_gamma))
# if there are invariant variables, then ncol(QX)<p, otherwise ncol(QX)=p
      
      if(vc.method=="SA"){
### SA var_u
        XPy<-crossprod(X_tilde, Py)
        v_b<-y_tilde-X_tilde%*%invXPX%*%XPy
        v_b_bar<-tapply(v_b, cls, mean)
        vPv<-sum(ni*v_b_bar*v_b_bar) 
        var_u<-max( (vPv-(m-p)*var_e)/(N-sum(diag(invXPX%*%sum_xixi))), 0)

### FC var_u 
      }else if(vc.method=="FC"){
        v<-lm(y_tilde~X_tilde-1)$residuals
        vv<-v%*%v
        var_u<-ifelse(p==1L, 
                     (vv-(N-p)*var_e)/(N-sum_nixx/XX),
                     (vv-(N-p)*var_e)/(N-sum(diag(solve(XX)%*%sum_nixx))) )
        var_u<-max(var_u, 0)
      }
    }
    
    if(vc.method=="FChetero"){
      alpha<-plmm$alpha
      ay<-y_tilde/alpha
      aX<-X_tilde/alpha
      v_bar<-lm(ay~aX-1)$residuals
      vv<-v_bar%*%v_bar 
      
      eQe_eta_Cfc<-FC_hetero(N=N, p=p, m=m, plmm=plmm, invarIND=invarIND)
      eQe<-eQe_eta_Cfc$eQe
      eta_Cfc<-eQe_eta_Cfc$eta
      
      var_e<-eQe/(N-m-ncolQX-df_gamma)
      var_u<-max( (vv-(N-p)*var_e)/(eta_Cfc), 0)
      
    }else if(vc.method=="SAhetero"){
      eQe_vPv_eta_Csa<-SA_hetero(N=N, p=p, m=m, plmm=plmm, invarIND=invarIND)
      eQe<-eQe_vPv_eta_Csa$eQe
      vPv<-eQe_vPv_eta_Csa$vPv
      eta_Csa<-eQe_vPv_eta_Csa$eta_Csa
    
      var_e<-eQe/(N-m-ncolQX-df_gamma)
      var_u<-max( (vPv-(m-p)*var_e)/eta_Csa, 0)
    }

    VC_new<-c(var_u=var_u, var_e=var_e)
    
### check the convergence
    if(r>1){# to prevent convergence without iteration
      VCdiff<-abs(VC_new-VC)/VC
      VC<-VC_new
      if(max(VCdiff, na.rm=T) < epsilon){break}
    }
    
### beta and gamma estimation ###
    beta<-est_beta(VC=VC_new, vc.method=vc.method, p=p, plmm)
#    yXb=plmm$y-plmm$X%*%beta
    yXb<-plmm$y-X_tilde%*%beta #X_tilde = plmm$X

    if(nonpar.bws=="h.select"){
      h<-h.select(y=yXb, x=T_mat, method="cv", poly.index=poly.index, ...)
      
    }else if(nonpar.bws=="hcv"){
      h<-ifelse(d==1, hcv(y=yXb, x=as.vector(T_mat), poly.index=poly.index, ...), hcv(y=yXb, x=T_mat, poly.index=poly.index, ...) )
      
    }else if(nonpar.bws=="GCV"){
      sd_T<-apply(T_mat,2,sd)      
      h_start<-ifelse(d==1, sd_T/2, 0.5)         
      log_h<-optimize(f=GCV, interval=log(c(h_start/8, h_start*4)), yXb=yXb, d=d, N=N, poly.index=poly.index, plmm=plmm, sd_T=sd_T)
      h<-exp(log_h$minimum)
      if(d==2){h<-h*sd_T}      

    }else if(nonpar.bws=="GCV.c"){
      sd_T<-apply(T_mat,2,sd)
      h_start<-ifelse(d==1, sd_T/2, 0.5)
      
      if(vc.method == "FC" || vc.method=="SA"){
        plmm$R.i<-lapply(plmm$clsLev, .R.i, VC=VC_new, plmm=plmm)
      }else{## FC/SAhetero
        plmm$sd_ij0<-split(sqrt((plmm$alpha^2)*var_e+var_u), cls)
        plmm$R.i_hetero<-lapply(plmm$clsLev, .R.i_hetero, var_u=var_u, plmm=plmm)  
      }      
      plmm$GCV_c<-NULL
      plmm$df_gamma<-NULL

      log_h<-optimize(f=.GCV_c, interval=log(c(h_start/8, h_start*4)), yXb=yXb, d=d, N=N, poly.index=poly.index, plmm=plmm, sd_T=sd_T)
      h<-exp(log_h$minimum)
      
      df_gamma<-plmm$df_gamma[which(plmm$GCV_c==min(plmm$GCV_c))][1]
# df_gamma may have more than one same values if some GCV_c values are the same       
      if(d==2){h<-h*sd_T}
    }  
#h=2*h #<Experiment> with h

###gamma update
    h<-h*scale.h
    gamma<-sm.regression(y=yXb, x=T_mat, h=h, eval.grid=F, eval.points=T_mat, poly.index=poly.index, display="none")$estimate 
    
### Update df_gamma
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
    }

    r<-r+1    
    VC_iter[r+1, ]<-VC_new
    beta_iter[r+1, ]<-beta
  }
#### end of while() ####
  
  VC_iter<-na.omit(VC_iter)[,]
  beta_iter<-na.omit(beta_iter)[,]

  return(list(beta=beta, VC_iter=VC_iter, beta_iter=beta_iter, gamma=gamma, h=h, df_gamma=df_gamma, iter=r))
}
