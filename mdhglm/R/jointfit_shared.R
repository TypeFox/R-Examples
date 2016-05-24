jointfit_shared <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
PhiFix=NULL,LamFix=NULL,structure="shared",mord=0,dord=1,Maxiter=200,convergence=1e-06,Init_Corr=NULL,EstimateCorrelations=TRUE) {
    N_model<-length(RespDist)    
    mc <- match.call()
    model<-NULL
    ## initial values for dispersion ##
    if (is.null(Init_Corr)) ww<-rep(1,N_model)
    else ww<-Init_Corr
    phi<-rep(1,N_model)
    if (!is.null(PhiFix)) phi<-PhiFix
    lambda_value<-0.5
    if (!is.null(LamFix)) lambda_value<-LamFix
    #########################
    for (iii in 1:N_model) {
       res1<-MakeModel(RespDist=RespDist[iii],DataMain=DataMain[[iii]],MeanModel=MeanModel[[iii]])
       if (iii==1) { 
          yy<-res1[[1]]
          xx<-res1[[2]]
          zz<-ww[iii]*res1[[3]]
##        namesXX<-res1[[4]]
          namesYY<-res1[[5]]
          nn<-res1[[6]]
          pp<-res1[[7]]
          qq<-res1[[8]]
          RespLink<-MeanModel[[iii]][2][[1]]
          RandDist<-MeanModel[[iii]][4][[1]]
       } else {
          yy<-rbind(yy,res1[[1]])
          xx<-dbind(xx,res1[[2]])
          zz<-rbind(zz,ww[iii]*res1[[3]])
##        namesXX<-cbind(namesXX,res1[[4]])
          namesYY<-cbind(namesYY,res1[[5]])
          nn<-cbind(nn,res1[[6]])
          pp<-cbind(pp,res1[[7]])
          qq<-cbind(qq,res1[[8]])
          RespLink<-cbind(RespLink,MeanModel[[iii]][2][[1]])
          RandDist<-cbind(RandDist,MeanModel[[iii]][4][[1]])
       }
    }
    cum_n<-cumsum(c(0,nn))
    cum_q<-cumsum(c(0,qq))
    cum_p<-cumsum(c(0,pp))
    ## initial values for beta ##
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       if (RespDist[iii]=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink[iii]))
       if (RespDist[iii]=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink[iii]))
       if (RespDist[iii]=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink[iii]))
       if (RespDist[iii]=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink[iii]))
       temp<-matrix(0,pp[iii],1)
       temp[1:pp[iii],1]<-c(resglm$coefficients)[1:pp[iii]]
       if (iii==1) {
            beta_init<-temp
            beta_h<-temp
       }
       else {
            beta_init<-dbind(beta_init,temp)
            beta_h<-rbind(beta_h,temp)
       }
    }
    ## Estimation of random effects
       beta_mu<-beta_init
       v_h<-matrix(0,qq[1],1)
       xbeta1<-matrix(0,cum_n[N_model+1],1)
       mu<-matrix(0,cum_n[N_model+1],1)
       detadmu<-matrix(0,cum_n[N_model+1],1)
       Vmu<-matrix(0,cum_n[N_model+1],1)
       off <- matrix(0,cum_n[N_model+1],1)
       disp_est <- matrix(0,cum_n[N_model+1],1)
       lambda<-matrix(lambda_value,qq[1],1)
       dhdv<-matrix(0,cum_n[N_model+1],1)
   iter_v<-5
   convergence3<-100
   iteration<-1
  while (convergence3>convergence && iteration<=Maxiter ) {
   for (kkk in 1:iter_v) {
       for (iii in 1:N_model) {
          temp1<-cum_p[iii]+1
          temp2<-cum_p[iii+1]
          beta_mu[temp1:temp2,iii]<-beta_h[temp1:temp2,1]
       }
       xbeta<-xx %*% beta_mu
       for (iii in 1:N_model) {
          temp1<-cum_n[iii]+1
          temp2<-cum_n[iii+1]
          xbeta1[temp1:temp2,1]<-xbeta[temp1:temp2,iii]
          disp_est[temp1:temp2,1] <- phi[iii]
       }
       eta_mu <- off+xbeta1 + zz %*% v_h
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
    if (RespLink[iii]=="identity") {
        mu[temp1:temp2,1] <- eta_mu[temp1:temp2,1]
        detadmu[temp1:temp2,1] <- (abs(mu[temp1:temp2,1])+1)/(abs(mu[temp1:temp2,1])+1)
    }
    if (RespLink[iii]=="log") {
        mu[temp1:temp2,1] <- exp(eta_mu[temp1:temp2,1])
        detadmu[temp1:temp2,1] <- 1/mu[temp1:temp2,1]
    }
    if (RespLink[iii]=="logit") {
        mu[temp1:temp2,1] <- 1/(1+exp(-eta_mu[temp1:temp2,1]))
        detadmu[temp1:temp2,1] <- 1/(mu[temp1:temp2,1]*(1-mu[temp1:temp2,1]))
    }
    if (RespLink[iii]=="probit") {
        mu[temp1:temp2,1] <- pnorm(eta_mu[temp1:temp2,1])
        detadmu[temp1:temp2,1] <- 1/dnorm(eta_mu[temp1:temp2,1])
    }
    if (RespLink[iii]=="cloglog") {
        mu[temp1:temp2,1] <- 1-exp(-exp(eta_mu[temp1:temp2,1]))
        detadmu[temp1:temp2,1] <- 1/exp(-exp(eta_mu[temp1:temp2,1]))
    }
    if (RespDist[iii]=="gaussian") Vmu[temp1:temp2,1]<-(abs(mu[temp1:temp2,1])+1)/(abs(mu[temp1:temp2,1])+1)
    if (RespDist[iii]=="poisson") Vmu[temp1:temp2,1]<-mu[temp1:temp2,1]
    if (RespDist[iii]=="binomial") Vmu[temp1:temp2,1]<-mu[temp1:temp2,1]*(1-mu[temp1:temp2,1])
    if (RespDist[iii]=="gamma") Vmu[temp1:temp2,1]<-mu[temp1:temp2,1]^2
    }
       dmudeta<-1/detadmu
       temp4<-dmudeta^2 /(disp_est*Vmu)
       vector_w1<-temp4
       W1<-diag(as.vector(temp4))
       z1<-eta_mu+(yy-mu)*detadmu-off
       W2<-diag(1/as.vector(lambda))
       dhdv<-t(zz)%*%W1%*%(detadmu*(yy-mu))-W2%*%v_h
       d2hdv2<--t(zz)%*%W1%*%zz-W2
       v_h_old<-v_h
       v_h<-v_h-solve(d2hdv2)%*%dhdv
       c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
 ##      print(c_v_h)
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    n<-cum_n[N_model+1]
    a<-matrix(0,n,1)
   if (mord==1) {
    III<-diag(rep(1,qq[1]))
    T<-t(cbind(t(zz),III))
    Null1<-matrix(0,n,qq[1])
    Null2<-matrix(0,qq[1],n)
    W<-matrix(0,n+qq[1],n+qq[1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qq[1])),]<-cbind(Null2,W2)   
    P<-T%*%solve(t(T)%*%W%*%T)%*%t(T)%*%W
    K1<--zz%*%solve(t(T)%*%W%*%T)%*%t(zz)
    K2<--solve(t(T)%*%W%*%T)
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    for (i in 1:n){
        d1[i]<-P[i,i]*detadmu[i]
        d2[i]<-0
        for (qqq in 1:n){
            d2[i]<-d2[i]+P[qqq,qqq]*K1[qqq,i]
        }
        if (RandDist[1]=="gaussian") d3[i]<-0
    }
    d<-d1+d2+d3
    s<-d*dmudeta/2
    a<-(solve(W1)+zz%*%solve(W2)%*%t(zz))%*%W1%*%(s*detadmu)
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    Sig<- zz %*% solve(W2) %*% t(zz) +solve(W1)
    invSig<-solve(Sig)
    beta_h<-solve(t(xx)%*%invSig%*%xx)%*%(t(xx)%*%invSig%*%(z1-a))
    se_beta<-sqrt(diag(solve(t(xx)%*%invSig%*%xx)))
    c_beta_h<-sum(abs(as.vector(beta_h_old)-as.vector(beta_h)))
##    print(c_beta_h)
    }
    deviance_residual<-matrix(0,cum_n[N_model+1],1)
    old_phi<-phi
    OO1<-matrix(0,qq[1],cum_p[N_model+1])
    Null1<-matrix(0,cum_n[N_model+1],qq[1])
    Null2<-matrix(0,qq[1],cum_n[N_model+1])
    III<-diag(rep(1,qq[1]))
    TT<-rbind(cbind(xx,zz),cbind(OO1,III))
    WW<-matrix(0,cum_n[N_model+1]+qq[1],cum_n[N_model+1]+qq[1])
    WW[c(1:cum_n[N_model+1]),]<-cbind(W1,Null1)
    WW[c((cum_n[N_model+1]+1):(cum_n[N_model+1]+qq[1])),]<-cbind(Null2,W2)   
    PP<-TT%*%solve(t(TT)%*%WW%*%TT)%*%t(TT)%*%WW
    se_log_phi<-rep(0,N_model)
    log_phi<-rep(0,N_model)
    log_phi<-log(phi)
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       temp5<-cum_q[iii]+1
       temp6<-cum_q[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       z<-zz[temp1:temp2,1:qq[iii]]
       n<-nn[iii]
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (RespDist[iii]=="gaussian") deviance_residual[temp1:temp2,1]<-(y-mu[temp1:temp2,1])^2
    if (RespDist[iii]=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual[temp1:temp2,1]<-2*y_zero*mu[temp1:temp2,1]+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu[temp1:temp2,1])-(y+0.00001-mu[temp1:temp2,1]))
    }
    if (RespDist[iii]=="binomial") deviance_residual[temp1:temp2,1]<-(1*(y==0))*2*log(1/(1-mu[temp1:temp2,1]))+(1*(y==1))*2*log(1/mu[temp1:temp2,1])
    if (RespDist[iii]=="gamma") deviance_residual[temp1:temp2,1]<-2*(-log(y/mu[temp1:temp2,1])+(y-mu[temp1:temp2,1])/mu[temp1:temp2,1])
    leverage<-rep(0,n)
    for (kk in 1:n) leverage[kk]<-PP[temp1+kk-1,temp1+kk-1]
    resp_disp<-deviance_residual[temp1:temp2]/(1-leverage)
    RespLink_disp<-"log"
    weight_disp<-(1-leverage)/2
##############################################################
######### GLM fit for phi #####################
##############################################################
   if (is.null(PhiFix)) {
      if (RespDist[iii]=="gaussian" || RespDist[iii]=="gamma") {
       x_disp<-matrix(1,n,1)
       resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weights=weight_disp)
       inv_disp<-1/resglm_disp$fitted.values
       disp_est[temp1:temp2,1]<-1/inv_disp
       phi[iii]<-disp_est[temp1,1]
       res2_1<-summary(resglm_disp,dispersion=2)
       log_phi[iii]<-log(phi[iii])
       se_log_phi[iii]<-res2_1$coefficients[2]
      }
    }
    }
##############################################################
######### GLM fit for lambda            #####################
##############################################################
    psi<-matrix(0,qq[1],1)
    old_lambda_est<-lambda
         if (RandDist[1]=="gaussian") {
              psi<-psi+0
              u_h<-v_h
              resp_lambda<-u_h^2
          }
          if (RandDist[1]=="gamma") {
              psi<-psi+1
              u_h<-exp(v_h)
              resp_lambda<-2*(1-log(u_h)-(1-u_h))
          }    
          if (RandDist[1]=="beta") {
              psi<-psi+0.5
              u_h<-1/(1+exp(-v_h))
              resp_lambda<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
          }
          if (RandDist[1]=="inverse-gamma") {
              psi<-psi+1
              u_h<-exp(v_h)
              resp_lambda<-2*(log(u_h)+(1-u_h)/u_h)
              resp_lambda<-(resp_lambda>0)*resp_lambda+(resp_lambda<=0)*0.0001
          }
    leverage_lambda<-rep(0,qq[1])
    for (kk in 1:qq[1]) leverage_lambda[kk]<-PP[cum_n[N_model]+kk,cum_n[N_model]+kk]
    resp_lambda<-resp_lambda/(1-leverage_lambda)
    weight_lambda<-(1-leverage_lambda)/2
    log_lambda<-log(lambda_value)
    se_log_lambda<-0
  if (is.null(LamFix)) {
       x_lambda<-matrix(1,qq[1],1)
       RespLink_lambda="log"
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda),weights=weight_lambda)
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       res3_1<-summary(resglm_lambda,dispersion=2)
       log_lambda<-res3_1$coefficients[1]
       se_log_lambda<-res3_1$coefficients[2]
       convergence2<-sum(abs(lambda_est[1]-old_lambda_est[1]))
    } else convergence2<-0
##############################################################
######### Estimation of shared parameters   ##################
##############################################################
    old_ww<-ww
    if(N_model>1) {
       OO1<-matrix(0,qq[1],cum_p[N_model+1])
       Null1<-matrix(0,cum_n[N_model+1],qq[1])
       Null2<-matrix(0,qq[1],cum_n[N_model+1])
       III<-diag(rep(1,qq[1]))
       TT<-rbind(cbind(xx,zz),cbind(OO1,III))
       WW<-matrix(0,cum_n[N_model+1]+qq[1],cum_n[N_model+1]+qq[1])
       WW[c(1:cum_n[N_model+1]),]<-cbind(W1,Null1)
       WW[c((cum_n[N_model+1]+1):(cum_n[N_model+1]+qq[1])),]<-cbind(Null2,W2)
       solve_TT<-solve(t(TT)%*%WW%*%TT)
       se_ww<-rep(0,N_model)  
       for (iii in 2:N_model) {
            temp1<-cum_n[iii]+1
            temp2<-cum_n[iii+1]
            temp3<-cum_p[iii]+1
            temp4<-cum_p[iii+1]
            temp5<-cum_q[iii]+1
            temp6<-cum_q[iii+1]
            y<-yy[temp1:temp2,1]
            x<-xx[temp1:temp2,temp3:temp4]
            z<-zz[temp1:temp2,1:qq[iii]]
            W_temp<-W1[temp1:temp2,temp1:temp2]
            dhdw<-t(v_h)%*%t(z/ww[iii])%*%W_temp%*%(detadmu[temp1:temp2,1]*(y-mu[temp1:temp2,1]))
##            print(dhdw)
            d2hdw2<- -t(v_h)%*%t(z/ww[iii])%*%W_temp%*%(z/ww[iii])%*%v_h
            zz1<-0*zz
            zz1[temp1:temp2,1:qq[iii]]<-zz[temp1:temp2,1:qq[iii]]/ww[iii]
            dTTdw<-rbind(cbind(0*xx,zz1),cbind(0*OO1,0*III))
            TT11<-t(dTTdw)%*%WW%*%TT
            TT12<-t(TT)%*%WW%*%dTTdw
            dhdw<-dhdw-0.5*sum(diag(solve_TT%*%(TT11+TT12)))
##            print(dhdw)
##            print(d2hdw2)
            if (EstimateCorrelations==TRUE) ww[iii]<-ww[iii]+dhdw/(-d2hdw2)
            if (EstimateCorrelations==TRUE) se_ww[iii]<-sqrt(1/(-d2hdw2))
            zz[temp1:temp2,1:qq[iii]]<-zz[temp1:temp2,1:qq[iii]]/old_ww[iii]*ww[iii]
        }
       convergence4<-sum(abs(ww-old_ww))
     } else convergence4<-0
       convergence1<-sum(abs(phi-old_phi))
       convergence3<-convergence1+convergence2+convergence4
##       print("phi")
##       print(old_phi)
##       print("lambda")
##       print(lambda[1])
##       print("shared parameter")
##       print(ww)
       print_err<-convergence3
       names(print_err) <- "convergence : "
       print_i<-iteration
       names(print_i) <- "iteration : "
##       print(print_i)
##       print(print_err)
       iteration<-iteration+1
    } ## for loop
###############################################################
############# likelihood estimates ############################
###############################################################
    pi<-3.14159265359
    d2hdv2<--t(zz)%*%W1%*%zz-W2
    H<-t(zz)%*%W1%*%zz+W2
    X<-xx
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%zz%*%III)),cbind((t(III)%*%t(zz)%*%W1%*%X),H))
    hlikeli<-0
    pvh<-0
    pbvh<-0
  for (iii in 1:N_model) {
    temp1<-cum_n[iii]+1
    temp2<-cum_n[iii+1]
    y<-yy[temp1:temp2,1]
    mu1<-mu[temp1:temp2,1]
    if (RespDist[iii]=="gaussian") hlikeli<-hlikeli+sum(-0.5*(y-mu1)*(y-mu1)/phi[iii])-0.5*nn[iii]*log(2*pi*phi[iii])
    if (RespDist[iii]=="poisson") hlikeli<-hlikeli+sum(y*log(mu1)-mu1-lgamma(y+1))
    if (RespDist[iii]=="binomial") hlikeli<-hlikeli+sum(y*log(mu1)+(1-y)*log(1-mu1))
    if (RespDist[iii]=="gamma") hlikeli<-hlikeli+sum(log(y)/phi-log(y)-y/(phi[iii]*mu1)-log(phi[iii])/phi[iii]-log(mu1)/phi[iii]-lgamma(1/phi[iii]))
   }
    AA<-rbind(cbind((t(xx)%*%W1%*%xx),(t(xx)%*%W1%*%zz)),cbind((t(zz)%*%W1%*%xx),(-1*d2hdv2)))
    BB<-rbind(cbind((t(xx)%*%W1%*%xx),(t(xx)%*%W1%*%zz)),cbind((t(zz)%*%W1%*%xx),(t(zz)%*%W1%*%zz)))
    pd<- sum(diag(solve(AA) %*% BB))    
    caic<--2*hlikeli+2*pd
    cc1<-svd(W2)
    logdet1<-sum(log(abs(1/cc1$d)))
    hlikeli<-hlikeli-0.5*t(v_h)%*%W2%*%v_h-0.5*logdet1-0.5*log(2*pi*nrow(W2))
    cc1<-svd(-d2hdv2)
    logdet1<-sum(log(abs(cc1$d)))
    pvh<-hlikeli-0.5*logdet1+0.5*log(2*pi*nrow(d2hdv2))
    cc1<-svd(A)
    logdet1<-sum(log(abs(cc1$d)))
    pbvh<-hlikeli-0.5*logdet1+0.5*log(2*pi*nrow(A))
    m2h<--2*hlikeli
    m2pvh<--2*pvh
    m2pbvh<--2*pbvh       
############################################################## 
    res<-list(Mu=mu,V_h=v_h,Beta_h=beta_h,SE_Beta=se_beta,Log_Phi=log_phi,SE_Log_Phi=se_log_phi,Log_Lambda=log_lambda,SE_Log_Lambda=se_log_lambda,
              Shared=ww,SE_Shared=se_ww,M2h=m2h,M2pvh=m2pvh,M2pbvh=m2pbvh,CAIC=caic)
    return(res)
}
