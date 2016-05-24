dhglmfit_run <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
PhiFix=NULL,LamFix=NULL,mord=0,dord=1,REML=TRUE,Maxiter=200,convergence=1e-02,Iter_mean=3) {
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
    x_lambda <- matrix(fr_lambda$X)
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    if (DispersionModel[3][[1]]=="constant") DispersionModel[3][[1]]<-phi~1
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if (RespDist=="poisson" || RespDist=="binomial") PhiFix<-1
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    beta_mu<-matrix(0,p,1)
    beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
    RandDist2<-rep(0,nrand)
    RandDist1<-MeanModel[4][[1]]
    check<-0
    length3<-length(RandDist1)
   if (length3>1) {
    if(nrand>1) {
    for (i in 1:nrand) {
       if (RandDist1[i]=="gaussian") RandDist2[i]<-1
       if (RandDist1[i]=="gamma") RandDist2[i]<-2
       if (RandDist1[i]=="inverse-gamma") RandDist2[i]<-3
       if (RandDist1[i]=="beta") RandDist2[i]<-4
       if (i>1) check<-check+abs(RandDist2[i]-RandDist2[i-1])
    }
    }
    }
    if (q[1]>0) {
       qcum <- cumsum(c(0, q))
       v_h<-matrix(0,qcum[nrand+1],1)
       se_v_h<-matrix(0,qcum[nrand+1],1)
       u_h<-matrix(1,qcum[nrand+1],1)
       if (nrand>1) {
          RandDist1<-MeanModel[4][[1]]
          RandDist<-RandDist1[1]
       } else RandDist<-MeanModel[4][[1]]
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          RandDist1<-MeanModel[4][[1]]
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      }        
       oq<-matrix(1,qcum[nrand+1],1)
       if (is.null(LamFix)) temp6<-0.5
       else temp6<-LamFix
       lambda<-matrix(temp6,qcum[nrand+1],1)
       old_lambda_est<-lambda
       alpha_h <- rep(temp6, nrand)
       for (i in 1:nrand) {
          index1<-qcum[i]+1
          lambda[index1:qcum[i+1]]<-alpha_h[i]
       } 
    }
    if (q_disp[1]>0) {
       qcum_disp <- cumsum(c(0, q_disp))
       v_h_disp<-matrix(0,qcum_disp[nrand+1],1)
       RandDist_disp<-DispersionModel[4][[1]]
       if (RandDist=="gaussian") u_h_disp <- v_h_disp
       if (RandDist=="gamma") u_h_disp <-exp(v_h_disp)
       if (RandDist=="inverse-gamma") u_h_disp <-exp(v_h_disp)
       oq_disp<-matrix(1,qcum_disp[nrand+1],1)
       temp7<-exp(-3.40)
       lambda_disp<-matrix(temp7,qcum_disp[nrand+1],1)
       alpha_h_disp <- rep(temp7, nrand_disp)
    }
    convergence<-1e-03
while (convergence3>convergence && max_iter<=Maxiter ) {
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (q[1]==0) {
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weights=matrix(inv_disp),offset=Offset)
       beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
       eta_mu <- off + x %*% beta_mu 
    } 
##############################################################
######### HGLM estimates for mu          #####################
##############################################################
    if (q[1]==0) Iter_mean<-1
    if (q[1]>0) Iter_mean<-Iter_mean
  for (j in 1:Iter_mean) {
    if (q[1]>0) eta_mu <- off + x %*% beta_mu + z %*% v_h 
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta_mu
        detadmu <- -1/mu^2
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta_mu+(y-mu)*detadmu-off
    beta_h<-beta_mu
##############################################################
############# random effect  #################################
##############################################################
  if(q[1]>0) {
    I<-diag(rep(1,qcum[nrand+1]))
    W2<-diag(1/as.vector(lambda))
    c_v_h<-1.0
    iter_v<-1
        for (iiiii in 1:iter_v) {
    eta <- off + x %*% beta_h + z %*% v_h 
    eta_mu <- eta
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
  if (check==0) {
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,y-mu)-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="gamma") {
            dhdv<-crossprod(z,y-mu)+1/lambda-exp(v_h)/lambda  
            temp5<-exp(v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda)+exp(-v_h)/lambda  
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    } else {

  dhdv<-matrix(0,qcum[nrand+1],1)
  d2hdv2<-matrix(0,qcum[nrand+1],qcum[nrand+1])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
    zzz1<-zzz[[i]]
    if (RespDist=="poisson") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)+1/lambda[temp11:temp12]-exp(v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(v_h[temp11:temp12])/lambda[temp11:temp12]
            W21<-diag(as.vector(temp5))
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="binomial") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gamma") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="inverse-gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda[temp11:temp12])+exp(-v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(-v_h[temp11:temp12])/lambda[temp11:temp12]
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
  }

   }
    v_h_old<-v_h
    v_h<-v_h+dhdv/diag(-d2hdv2)
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      } 
                   }
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    a<-matrix(0,n,1)
    s<-matrix(0,n,1)
    if (q[1]>0 && mord==1 && RespDist!="gamma" && RespDist!="gaussian") {
    T<-t(cbind(t(matrix(z,nrow(z),ncol(z))),matrix(I,nrow(I),ncol(I))))
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)    
    K2<--solve(crossprod(T,diag(W)*T))  
    P<-T%*%(-K2)%*%crossprod(T,W)  
    K1<--z%*%(-K2)%*%t(z)  
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    if (RespDist=="binomial") d1<-diag(P)[1:n]*detadmu*(1-2*mu/BinomialDen)
    else d1<-diag(P)[1:n]*detadmu
    if (RandDist=="gaussian") d3<-0
     one<-matrix(1,ncol(K1),1)
     if (RespDist=="binomial") d2<-K1%*%((diag(P)[1:n]*(1-2*mu/BinomialDen))*one)
     else d2<-K1%*%(diag(P)[1:n]*one)
    d<-as.vector(d1)+as.vector(d2)+as.vector(d3)
    s<-d*dmudeta/2
    a<-(diag(1/diag(W1))+z%*%(1/diag(W2)*t(z)))%*%(diag(W1)*(s*detadmu))  
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    Sig<- z %*% (1/diag(W2) * t(z)) + diag(1/diag(W1))  
    invSig<-solve(Sig)  
    solve_xsx<-solve(crossprod(x,invSig%*%x))
    beta_h<-solve_xsx%*%(crossprod(x,invSig%*%(z1-a)))  
    se_beta<-sqrt(diag(solve_xsx))  
    beta_mu<-beta_h
############################################################## 
  }
}
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (q[1]==0) {
        diag<-glm.diag(resglm)
        leverage<-diag$h
    }
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
    if (q[1]>0) {
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
       TT<-rbind(cbind(x,z),cbind(OO1,I))
       WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW[c(1:n),]<-cbind(W1,Null1)
       WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
       K2<-solve(crossprod(TT,diag(WW)*TT))
       PP<-TT%*%K2%*%t(diag(WW)*TT)   
       leverage<-rep(0,n)
       leverage[1:n]<-diag(PP)[1:n]
    }
    resp_disp<-deviance_residual/(1-leverage)
    resp_disp_zero<-(resp_disp>0)*1
    resp_disp<-resp_disp_zero*resp_disp+(1-resp_disp_zero)*0.001
    RespLink_disp<-DispersionModel[2][[1]]
    Offset_disp<-DispersionModel[5][[1]]
    weight_disp<-(1-leverage)/2

##############################################################
######### GLM fit for phi #####################
##############################################################
  if (is.null(PhiFix)) {
    if (q_disp[1]==0) {
      if (RespDist=="gaussian" || RespDist=="gamma") {
#       resglm_disp<-glm(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1,family=Gamma(link=RespLink_disp),weights=weight_disp,offset=Offset_disp)
       resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weights=weight_disp,offset=Offset_disp)
       inv_disp<-1/resglm_disp$fitted.values
       disp_est<-1/inv_disp
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
      }
    }
##############################################################
######### HGLM fit for phi #####################
##############################################################
    if (q_disp[1]>0) {
       model_number=4
       RandDist_disp<-DispersionModel[4][[1]]
       disp_rand<-FL_disp$Subject[[1]]
       DataMain1<-list(matrix(resp_disp),matrix(x_disp),matrix(disp_rand))
       reshglm_disp<-hglmfit_corr(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1+(1|disp_rand),DataMain=DataMain1,Offset=Offset_disp,RespDist="gamma",
                                RespLink=RespLink_disp,RandDist=RandDist_disp,Maxiter=2,Iter_mean=2)
       disp_est<-reshglm_disp[10][[1]]
       inv_disp<-1/reshglm_disp[10][[1]]
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
    }
  } else convergence1<-0
    if (q[1]>0) {
       z_dimension<-rep(0,nrand)
       for (i in 1:nrand) z_dimension[i]<-qcum[i+1]-qcum[i]
       psi<-matrix(0,qcum[nrand+1],1)
       resp_lambda<-matrix(0,qcum[nrand+1],1)
       leverage1<-rep(0,qcum[nrand+1])
     if (check==0) {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist=="gaussian") {
              psi<-psi+0
              temp17<-u_h^2
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="gamma") {
              psi<-psi+1
               temp17<-2*(-log(u_h)-(1-u_h))
 #             temp17<-2*(-log(u_h)-(1-u_h)+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
               temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1])
 #             temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1]+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }    
          if (RandDist=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h)+(1-u_h)/u_h)
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
       }
    } else {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist1[i]=="gaussian") {
              psi<-psi+0
              temp17<-u_h[temp16:qcum[i+1]]^2
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="gamma") {
              psi<-psi+1
              temp17<-2*(-log(u_h[temp16:qcum[i+1]])-(1-u_h[temp16:qcum[i+1]]))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }    
          if (RandDist1[i]=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h[temp16:qcum[i+1]])+(1-0.5)*log((1-0.5)/(1-u_h[temp16:qcum[i+1]])))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h[temp16:qcum[i+1]])+(1-u_h[temp16:qcum[i+1]])/u_h[temp16:qcum[i+1]])
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
       }
    }
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
##       TT<-rbind(cbind(x,z),cbind(OO1,I))
##       WW<-rbind(cbind(W1,Null1),cbind(Null2,W2))
       WW2<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW2[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2) 
##       solve_twt<-solve(crossprod(TT,diag(WW)*TT))
       solve_twt<-K2
##       PP<-TT%*%solve_twt%*%t(diag(WW)*TT)   
       TT1<-rbind(z,I)
       if (REML==FALSE) {
             solve_t1wt1<-solve(crossprod(TT1,diag(WW)*TT1))
             PP<-TT1%*%solve_t1wt1%*%t(diag(WW)*TT1)
       }
       solve_zw1zw2<-solve(crossprod(z,diag(W1)*z)+W2)
##       dvdlam<- solve_zw1zw2%*%(diag(W2)*diag(W2)*v_h)
       dvdlam1<- solve_zw1zw2%*%(diag(W2)*v_h)
       if (RespDist=="binomial") temp5<-(1-2*mu/BinomialDen)*mu/BinomialDen*(1-mu/BinomialDen)
       if (RespDist=="poisson") temp5<-mu
       if (RespDist=="gamma" || RespDist=="gaussian") temp5<-0*mu
       dWdmu<-diag(as.vector(temp5))
       dWdlam<-diag(as.vector((diag(dWdmu)*z)%*%dvdlam1))
       dWWdlam<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       dWWdlam[c(1:n),]<-cbind(dWdlam,Null1)
       HHstar<-crossprod(TT,diag(dWWdlam)*TT)
       AAA<-solve_twt%*%HHstar
       dWW1dlam<-cbind(dWdlam,Null1)
       if (REML==FALSE) {
             HHstar1<-crossprod(TT1,diag(dWWdlam)*TT1)
             AAA<-solve_t1wt1%*%HHstar1
       }
       aaaa<-rep(0,qcum[nrand+1])
       if (RespDist=="binomial") adj<-2*(crossprod(z,y-mu)-diag(W2)*v_h)*dvdlam1
       else adj<-0
       index1<-qcum[nrand+1]
       index2<-n+1
       index3<-n+qcum[nrand+1]
       index4<-p+1
       index5<-p+qcum[nrand+1]
       if (REML==TRUE) leverage1[1:index1]<-diag(PP)[index2:index3]
       else leverage1[1:index1]<-diag(PP)[index2:index3]
       if (REML==TRUE) aaaa[1:index1]<--diag(AAA)[index4:index5]-sum(diag(AAA)[1:p])/qcum[nrand+1]
       else aaaa[1:index1]<--diag(AAA)[1:index1]
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           pppp<-mu/BinomialDen
	   dh5dv5<-(1-14*pppp+36*pppp*pppp-24*pppp*pppp*pppp)*pppp*(1-pppp)
	   dh4dv4<-(1-6*pppp+6*pppp*pppp)*pppp*(1-pppp)
	   dh3dv3<-(1-2*pppp)*pppp*(1-pppp)
	   dh2dv2<-pppp*(1-pppp)
           d4hdv4dlam<--t(z)%*%diag(as.vector(dh5dv5))%*%z%*%diag(as.vector(dvdlam1))
           d3hdv3dlam<--t(z)%*%diag(as.vector(dh4dv4))%*%z%*%diag(as.vector(dvdlam1))
           d2hdv2dlam<-t(z)%*%diag(as.vector(dh3dv3))%*%z%*%diag(as.vector(dvdlam1))-W2
           d2hdv2<-t(z)%*%(diag(W1)*z)+W2
           d4hdv4<--t(z)%*%diag(as.vector(dh4dv4))%*%z
           d3hdv3<--t(z)%*%diag(as.vector(dh3dv3))%*%z
           bbb1<-3/24*d4hdv4dlam%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb2<--3/24*2*d4hdv4%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb3<-5/24*2*d3hdv3dlam%*%d3hdv3%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb4<--5/24*3*d3hdv3%*%d3hdv3%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb<-(bbb1+bbb2+bbb3+bbb4)
       }
       bbbb<-rep(0,qcum[nrand+1])
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           bbbb<- 1*diag(2*bbb)
       }
##       leverage1<-0
###       print(mean(u_h^2/(1-leverage1-aaaa)))
###       aaaa<-0
       if (REML==TRUE) resp_lambda<-(resp_lambda)/(1-leverage1-aaaa-bbbb)
       else resp_lambda<-resp_lambda/(1-leverage1-aaaa-bbbb)
       resp_lambda_neg<-1*(resp_lambda<0)
       resp_lambda<-(1-resp_lambda_neg)*resp_lambda+resp_lambda_neg*0.0001
       weight_lambda<-abs((1-leverage1-aaaa-bbbb)/2)
##       if (nrand==3 && check!=0) {
##           resp_lambda<-resp_lambda*(1-leverage)
##           weight_lambda<-weight_lambda/weight_lambda
##       }
    }
     maximum<-10
     if (nrand>=3) maximum<-5
##############################################################
######### GLM fit for lambda            #####################
##############################################################
 if (length1<=1) {
  if (is.null(LamFix)) {
    if (q[1]>0 && q_lambda[1]==0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
##       print(resp_lambda)
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       tttt<-sum(lambda_est/lambda_est)
##       convergence2<-sum(abs(lambda_est-old_lambda_est))/tttt
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence2<-0
  } else convergence2<-0
    convergence3<-convergence1+convergence2
    if (model_number==1) convergence3<-0
##    print_i<-max_iter
##    print_err<-convergence3
##    names(print_i) <- "iteration : "
##    print(print_i)
##    names(print_err) <- "convergence : "
##    print(print_err)
##    max_iter<-max_iter+1
## }
##############################################################
######### HGLM fit for lambda            #####################
##############################################################
    if (q[1]>0 && q_lambda[1]>0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       RespLink_lambda<-MeanModel[7][[1]]
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       RandDist_lambda<-MeanModel[9][[1]]
       RespLink_lambda<-MeanModel[7][[1]]
       x_lambda<-matrix(1,q_lambda[1],1)
       lambda_rand<-c(1:q_lambda[1])
       resp_lambda1<-resp_lambda[1:q_lambda[1]]
       resp_lambda<-resp_lambda1
       DataMain2<-list(resp_lambda,x_lambda,lambda_rand)
       model_number1<-1
       reshglm_lambda<-hglmfit_corr(resp_lambda~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[1:nnn]<-lambda_est1[1:nnn,1]
       lambda_est<-lambda
##       convergence21<-sum(abs(lambda_est-old_lambda_est))/nnn
       convergence21<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence21<-0
  }
   if(length1>1) {

  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      if (nrand>=3 && p_lambda>5) indicator1<-0
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       convergence21<-convergence21+sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weights=weight_lambda1)
       aaa<-summary(resglm_lambda)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
       convergence2<-convergence2+sum(abs(lambda_est-old_lambda_est))
     }
   }
    } ## length1
##    print(convergence1)
##    print(convergence2)
##    print(convergence21)
    convergence3<-convergence1+convergence2+convergence21
    print_i<-max_iter
    print_err<-convergence3
    names(print_i) <- "iteration : "
#    print(print_i)
    names(print_err) <- "convergence : "
#    print(print_err)
    max_iter<-max_iter+1
}
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/sqrt(1-leverage)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/(sqrt(1-leverage))
    md<-RespDist
    names(md)<-"Distribution of Main Response : "
    print(md)
    print("Estimates from the model(mu)")
    print(formulaMean)
    print(RespLink)
#    print(mean_residual)
#   print(reshglm_lambda)
    if (q[1]==0) {
        res1<-summary(resglm)
        beta_h<-beta_mu
        temp14<-p+1
        temp15<-2*p
        se_beta<-res1$coefficients[temp14:temp15]
        z_beta<-beta_h/se_beta
##        pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
    }
    if (length1<=1) {
    if (q[1]>0 && q_lambda[1]==0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates for logarithm of lambda=var(u_mu)")
        print(MeanModel[4][[1]])
        res3<-summary(resglm_lambda,dispersion=2)
        p_lambda<-nrand
        if (nrand>=3 && p_lambda>5) indicator1<-0
        lambda_h<-res3$coefficients[1:p_lambda]
        lambda_h[1]<-lambda_h[1]
#        lambda_h[2]<-lambda_h[2]*1.7
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
#        lambda_se[2]<-lambda_se[2]*sqrt(1.7)
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesRE
        print(lambda_coeff,4)
    }        
    if (q_lambda[1]>0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(formulaLambda)
        print(RandDist_lambda)
        res5<-reshglm_lambda
        temp9<-p_lambda+1
        temp10<-2*p_lambda
        beta_lambda<-res5[2][[1]]
        se_lambda<-res5[3][[1]]
        z_lambda<-beta_lambda/se_lambda
        res3<-summary(resglm_lambda,dispersion=2)
        if (nrand==1) {
           lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesX_lambda
        }
        if (nrand>1) {
           lambda_h<-res3$coefficients[1:nrand]
           temp11<-nrand+1
           temp12<-2*nrand
           lambda_se<-res3$coefficients[temp11:temp12]
           z_lambda<-lambda_h/lambda_se
           lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesRE
        }
        print(lambda_coeff,4)
        print("Estimates for logarithm of alpha=var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
    }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of tau=var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
    }
    if (length1>1) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(MeanModel[4][[1]])
  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      if (iiii==1) print("Estimates from the model(lambda1=var(u_mu))")
      if (iiii==2) print("Estimates from the model(lambda2=var(u_mu))")
      if (iiii==3) print("Estimates from the model(lambda3=var(u_mu))")
      if (iiii==4) print("Estimates from the model(lambda4=var(u_mu))")
      if (iiii==5) print("Estimates from the model(lambda5=var(u_mu))")
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      random_lambda<-findbars(formulaLambda)
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
              RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       res5<-reshglm_lambda
       temp9<-p_lambda+1
       temp10<-2*p_lambda
       beta_lambda<-res5[2][[1]]
       se_lambda<-res5[3][[1]]
       z_lambda<-beta_lambda/se_lambda
       lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
       colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
       rownames(lambda_coeff) <- namesX_lambda
       print(lambda_coeff,4)
        print("Estimates for logarithm of alpha=var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
        model_number1<-1
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weights=weight_lambda1)
       res3<-summary(resglm_lambda,dispersion=2)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
        lambda_h<-res3$coefficients[1:p_lambda]
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesX_lambda
        print(lambda_coeff,4)
     }
   }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of tau=var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
   }
#    v_h1<-corr_res[[7]]
    pi<-3.14159265359
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-lgamma(y+1))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu/BinomialDen)+(BinomialDen-y)*log(1-mu/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/disp_est-log(y)-y/(disp_est*mu)-log(disp_est)/disp_est-log(mu)/disp_est-lgamma(1/disp_est))
    hh<-hlikeli
    if (RespDist=="gaussian") ll_y<-sum(-0.5*(y-y)*(y-y)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") ll_y<-sum((y+0.00001)*log(y)-y-lgamma(y+1))
    if (RespDist=="binomial") ll_y<-sum(y*log((y+0.00001)/BinomialDen)+(BinomialDen-y)*log(1-(y-0.00001)/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") ll_y<-sum(log(y)/disp_est-log(y)-y/(disp_est*y)-log(disp_est)/disp_est-log(y)/disp_est-lgamma(1/disp_est))
    if (RespDist=="gaussian") deviance<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") deviance<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    if (RespDist=="gamma") deviance<-2*(-log(y/mu)+(y-mu)/mu)
    if (RespDist=="gaussian" || RespDist=="gamma") deviance<-deviance/disp_est
    print("========== Likelihood Function Values and Condition AIC ==========")
    if (model_number == 1 || model_number == 2) {
       ml<- -2*hlikeli
       d2hdx2<--t(x)%*%(diag(W1)*x)  
       rl<- ml+log(abs(det(-d2hdx2/(2*pi))))
       pd<- p
       caic<-ml+2*pd
       sd<- -2*hh + 2*ll_y
       df<-length(y)-pd
    }
    if (model_number >=3) {
     if (check==0) {
       if (RandDist=="gaussian") {
            cc1<-svd(W2)
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*logdet1  
       } 
##       if (RandDist=="gaussian") hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2))))  
       if (RandDist=="gamma") hv<-log(u_h)/lambda_est-u_h/lambda-log(lambda_est)/lambda_est-lgamma(1/lambda_est)
       if (RandDist=="inverse-gamma") {
           lambda_est1<-lambda_est/(1+lambda_est)
           alpha<-(1-lambda_est1)/lambda_est1
###           hv<-(v_h-log(u_h))/lambda_est1-(1+1/lambda_est1)*log(lambda_est1)-lgamma(1/lambda_est1)+log(lambda_est1)
           hv<-(alpha+1)*(log(alpha)-v_h)-alpha/u_h-lgamma(alpha+1)
       }
       if (RandDist=="beta") {
           lambda_est1<-2*lambda_est/(1-lambda_est)
           hv<-(0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1)
       }
     } else {
    hv<-matrix(0,1,1)
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
       if (RandDist1[i]=="gaussian") {
            cc1<-svd(W2[temp11:temp12,temp11:temp12])
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<-hv-0.5*t(v_h[temp11:temp12])%*%W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]-0.5*nrow(W2[temp11:temp12,temp11:temp12])*log(2*pi)-0.5*logdet1  
       } 
       if (RandDist1[i]=="gamma") {
            hv<-hv+sum(log(u_h[temp11:temp12])/lambda_est[temp11:temp12]-u_h[temp11:temp12]/lambda[temp11:temp12]-log(lambda_est[temp11:temp12])/lambda_est[temp11:temp12]-lgamma(1/lambda_est[temp11:temp12]))
       }
       if (RandDist1[i]=="inverse-gamma") {
           lambda_est1<-lambda_est[temp11:temp12]/(1+lambda_est[temp11:temp12])
           alpha<-(1-lambda_est1)/lambda_est1
           hv<-hv+sum((alpha+1)*(log(alpha)-v_h[temp11:temp12])-alpha/u_h[temp11:temp12]-lgamma(alpha+1))
       }
       if (RandDist1[i]=="beta") {
           lambda_est1<-2*lambda_est[temp11:temp12]/(1-lambda_est[temp11:temp12])
           hv<-hv+sum((0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1))
       }
   }   
     }
       if (model_number == 4) hv10<-reshglm_disp[[8]][1]
       else hv10<-0
       if (model_number1 == 1) hv20<-reshglm_lambda[[8]][1]
       else hv20<-0
        if (model_number == 4) hv11<-reshglm_disp[[8]][2]
       else hv11<-0
       if (model_number1 == 1) hv21<-reshglm_lambda[[8]][2]
       else hv21<-0
       if (model_number == 4) {
            hv12<-reshglm_disp[[8]][3]
       }
       else hv12<-0
       if (model_number1 == 1) hv22<-reshglm_lambda[[8]][3]
       else hv22<-0
        cc1<-svd((-d2hdv2)/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       ml<- -2*hlikeli-2*sum(hv)+logdet1 ##-log(2*pi*nrow(d2hdv2))
##       ml<- -2*hlikeli-2*sum(hv)+log(abs(det(-d2hdv2/(2*pi))))
       W1x<-diag(W1)*x
       W1z<-diag(W1)*z
       AA<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((-1*d2hdv2),nrow(d2hdv2),ncol(d2hdv2))))  
       BB<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((t(z)%*%W1z),nrow(t(z)%*%W1z),ncol(t(z)%*%W1z))))  
        cc1<-svd(AA/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       rl<--2*hlikeli-2*sum(hv)+logdet1 
##       rl<--2*hlikeli-2*sum(hv)+logdet1-log(2*pi*nrow(AA))
       pd<- sum(diag(solve(AA) %*% BB))  
       caic<- -2*hlikeli + 2*pd
       sd<- -2*hh + 2*ll_y
       sd<-sum(deviance)
       df<-length(y)-pd
    }
#    likeli_coeff<-rbind(matrix(ml),matrix(rl),matrix(caic),matrix(sd),matrix(df))
    likeli_coeff<-rbind(ml,rl,caic,sd,df)
    if (model_number == 1 || model_number == 2) {
       rownames(likeli_coeff)<-c("-2ML (-2 h)          : ","-2RL (-2 p_beta (h)) : ","cAIC                 : ", "SD                   : ","df                   : ")
    }
    if (model_number ==3 ) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu) (h))          : ","-2RL (-2 p_beta(mu),v(mu) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
    if (model_number == 4) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
    if (model_number<4 && model_number1 == 1) {
        rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
    if (model_number == 4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
##    df<-n-p
##    Devdf<-matrix(c(df,sum(deviance_residual)),nrow=2)
##    rownames(Devdf)<-c("DF :","Deviance :")
    print(likeli_coeff)
	res <- list(mean_residual=matrix(mean_residual),mu=matrix(mu),deviance_residual=matrix(deviance_residual),df=df, ml = ml, rl = rl, caic = caic, md = md, RespLink = RespLink, beta_coeff = beta_coeff, RespLink_disp = RespLink_disp, phi_coeff = phi_coeff,alpha_coeff=alpha_coeff, tau_coeff=tau_coeff, RandDist = RandDist, lambda_coeff = lambda_coeff, scaled_dv = sd, df = df)
    return(res)
}
