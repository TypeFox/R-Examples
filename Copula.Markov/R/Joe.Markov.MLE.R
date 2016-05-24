Joe.Markov.MLE <-
function(Y,k=3,D=1,plot=TRUE){
  n=length(Y)   ##sample size##
  G=function(y,mu,sigma){pnorm((y-mu)/sigma)}   #G function
  g=function(y,mu,sigma){dnorm((y-mu)/sigma)}   #g function

  ##############Log-likelihood function #################
  L_function=function(theta){
    mu=theta[1]
    sigma=theta[2]
    alpha=theta[3]
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    Z=log(alpha-1+A)+(alpha-1)*log(1-U_t_1)+(alpha-1)*log(1-U_t)+(1/alpha-2)*log(A)
    ZZ=log(g(Y[1:n],mu,sigma)/sigma)
    -(sum(Z)+sum(ZZ))/n
  }
  #res=nlm(L_function,c(mean(Y),sd(Y),1),hessian=TRUE)
  #est=res$estimate
  
  ##############dL/dmu############################
  dL_dmu=function(mu,sigma,alpha){
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    A_mu1=alpha*(1-U_t_1)^(alpha-1)*(1-(1-U_t)^alpha)*u_t_1/sigma
    A_mu2=alpha*(1-U_t)^(alpha-1)*(1-(1-U_t_1)^alpha)*u_t/sigma
    A_mu=A_mu1+A_mu2
    L1=(Y[1:n]-mu)/sigma^2
    L2=A_mu/(alpha-1+A)
    L3=(alpha-1)/sigma*( u_t_1/(1-U_t_1)+u_t/(1-U_t) )
    L4=(1/alpha-2)*A_mu/A
    return((sum(L1)+sum(L2)+sum(L3)+sum(L4))/n)
  }

  ############dL/dsigma###########################
  dL_dsigma=function(mu,sigma,alpha){
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    A_s1=alpha*(1-U_t_1)^(alpha-1)*(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)*u_t_1/sigma^2
    A_s2=alpha*(1-U_t)^(alpha-1)*(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)*u_t/sigma^2
    A_s=A_s1+A_s2
    L1=(Y[1:n]-mu)^2/sigma^3
    L2=A_s/(alpha-1+A)
    L3=(alpha-1)/sigma^2*( (Y[1:n-1]-mu)*u_t_1/(1-U_t_1)+(Y[2:n]-mu)*u_t/(1-U_t) )
    L4=(1/alpha-2)*A_s/A
    return(-1/sigma+(sum(L1)+sum(L2)+sum(L3)+sum(L4))/n)
  }
  
  ############dL/dalpha#############################
  dL_dalpha=function(mu,sigma,alpha){
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    A_a1=(1-U_t_1)^alpha*log(1-U_t_1)+(1-U_t)^alpha*log(1-U_t)
    A_a2=(1-U_t_1)^alpha*(1-U_t)^alpha*log((1-U_t_1)*(1-U_t))
    A_a=A_a1-A_a2
    L1=log(1-U_t_1)+log(1-U_t)
    L2=(1+A_a)/(alpha-1+A)
    L3=-log(A)/alpha^2
    L4=(1/alpha-2)*A_a/A
    return((sum(L1)+sum(L2)+sum(L3)+sum(L4))/n)
  }
  
  ############F function##############################
  F=function(mu,sigma,alpha){
    c(dL_dmu(mu,sigma,alpha),dL_dsigma(mu,sigma,alpha),dL_dalpha(mu,sigma,alpha))
  }
  #F(est[1],est[2],est[3])
  
  ############d^2L/dmu^2 (in progress) ###########################
  d2L_dmu2=function(mu,sigma,alpha){
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    h_t_1=u_t_1/(1-U_t_1);h_t=u_t/(1-U_t)
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    
    A_mu1=alpha*(1-U_t_1)^(alpha-1)*(1-(1-U_t)^alpha)*u_t_1/sigma
    A_mu2=alpha*(1-U_t)^(alpha-1)*(1-(1-U_t_1)^alpha)*u_t/sigma
    A_mu=A_mu1+A_mu2
    
    B_mm1=(alpha-1)/(1-U_t_1)*u_t_1*(1-(1-U_t)^alpha)
    B_mm2=-alpha*(1-U_t)^(alpha-1)*u_t
    B_mm3=(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)/sigma
    A_mm1=alpha*(1-U_t_1)^(alpha-1)*(B_mm1+B_mm2+B_mm3)*u_t_1/sigma^2
    C_mm1=(alpha-1)/(1-U_t)*u_t*(1-(1-U_t_1)^alpha)
    C_mm2=-alpha*(1-U_t_1)^(alpha-1)*u_t_1
    C_mm3=(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)/sigma
    A_mm2=alpha*(1-U_t)^(alpha-1)*(C_mm1+C_mm2+C_mm3)*u_t/sigma^2
    A_mm=A_mm1+A_mm2
    
    L1=A_mm/(alpha-1+A)-( A_mu/(alpha-1+A) )^2
    L2=(alpha-1)/sigma^2*( (Y[1:n-1]-mu)*h_t_1/sigma-h_t_1^2 )
    L3=(alpha-1)/sigma^2*( (Y[2:n]-mu)*h_t/sigma-h_t^2 )
    L4=(1/alpha-2)*( A_mm/A-(A_mu/A)^2 )
    return((sum(L1)+sum(L2)+sum(L3)+sum(L4))/n-1/sigma^2)
  }
  
  ############d^2L/dsigma^2 (in progress)#####################
  d2L_dsigma2=function(mu,sigma,alpha){
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    h_t_1=u_t_1/(1-U_t_1);h_t=u_t/(1-U_t)
    
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    
    A_mu1=alpha*(1-U_t_1)^(alpha-1)*(1-(1-U_t)^alpha)*u_t_1/sigma
    A_mu2=alpha*(1-U_t)^(alpha-1)*(1-(1-U_t_1)^alpha)*u_t/sigma
    A_mu=A_mu1+A_mu2
       
    A_s1=alpha*(1-U_t_1)^(alpha-1)*(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)*u_t_1/sigma^2
    A_s2=alpha*(1-U_t)^(alpha-1)*(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)*u_t/sigma^2
    A_s=A_s1+A_s2
    
    B_ms1=(alpha-1)/(1-U_t_1)*u_t_1*(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)
    B_ms2=-alpha*(1-U_t)^(alpha-1)*u_t*(Y[1:n-1]-mu)-1+(1-U_t)^alpha
    B_ms3=(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)^2/sigma
    A_ms1=alpha*(1-U_t_1)^(alpha-1)*(B_ms1+B_ms2+B_ms3)*u_t_1/sigma^3
    C_ms1=(alpha-1)/(1-U_t)*u_t*(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)
    C_ms2=-alpha*(1-U_t_1)^(alpha-1)*u_t_1*(Y[2:n]-mu)-1+(1-U_t_1)^alpha
    C_ms3=(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)^2/sigma
    A_ms2=alpha*(1-U_t)^(alpha-1)*(C_ms1+C_ms2+C_ms3)*u_t/sigma^3
    A_ms=A_ms1+A_ms2
        
    L1=-2*(Y[1:n]-mu)/sigma^3
    L2=A_ms/(alpha-1+A)-A_mu*A_s/(alpha-1+A)^2
    L3=-(alpha-1)*( h_t_1/sigma^2+h_t/sigma^2 )
    L4=(alpha-1)*(  (Y[1:n-1]-mu)/sigma^3*( (Y[1:n-1]-mu)/sigma*h_t_1-u_t_1^2  ) )
    L5=(alpha-1)*(  (Y[2:n]-mu)/sigma^3*( (Y[2:n]-mu)/sigma*h_t-u_t^2  ) )
    #L6=(1/alpha-2)*( A_ss/A-(A_s/A)^2 )
  
    return(1/sigma^2+(  sum(L1)+sum(L2)+sum(L3)+sum(L4)+sum(L5)  )/n)
  }

  d2L_dsigma2=function(mu,sigma,alpha){
    h=0.0000001
    (dL_dsigma(mu,sigma+h,alpha)-dL_dsigma(mu,sigma,alpha))/h
  }
  
  #d2L_dsigma2(est[1],est[2],est[3])
  #-res$hessian[2,2]
  #h=0.0000001
  #-(L_function(est+2*c(0,h,0))-2*L_function(est+c(0,h,0))+L_function(est))/h^2
  
  ############d^2L/dalpha^2#################################
  d2L_dalpha2=function(mu,sigma,alpha){
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    A_a1=(1-U_t_1)^alpha*log(1-U_t_1)+(1-U_t)^alpha*log(1-U_t)
    A_a2=(1-U_t_1)^alpha*(1-U_t)^alpha*log((1-U_t_1)*(1-U_t))
    A_a=A_a1-A_a2
    A_aa1=(1-U_t_1)^alpha*(log(1-U_t_1))^2+(1-U_t)^alpha*(log(1-U_t))^2
    A_aa2=(1-U_t_1)^alpha*(1-U_t)^alpha*(log((1-U_t_1)*(1-U_t)))^2
    A_aa=A_aa1-A_aa2
    L1=A_aa/(alpha-1+A)-((1+A_a)/(alpha-1+A))^2
    L2=2*log(A)/alpha^3-2*A_a/A/alpha^2
    L3=(1/alpha-2)*(A_aa/A-(A_a/A)^2)
    return((sum(L1)+sum(L2)+sum(L3))/n)
  }

  ############d^2L/dmudsigma (in progress) #################################
  d2L_dmudsigma=function(mu,sigma,alpha){
 
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    h_t_1=u_t_1/(1-U_t_1);h_t=u_t/(1-U_t)
    
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    A_s1=alpha*(1-U_t_1)^(alpha-1)*(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)*u_t_1/sigma^2
    A_s2=alpha*(1-U_t)^(alpha-1)*(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)*u_t/sigma^2
    A_s=A_s1+A_s2
     
    B_ss1=(alpha-1)/(1-U_t_1)*u_t_1*(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)^2
    B_ss2=-alpha*(1-U_t)^(alpha-1)*u_t*(Y[2:n]-mu)*(Y[1:n-1]-mu)
    B_ss3=(1-(1-U_t)^alpha)*(Y[1:n-1]-mu)^3/sigma
    A_ss1=alpha*(1-U_t_1)^(alpha-1)*(B_ss1+B_ss2+B_ss3)*u_t_1/sigma^4
    C_ss1=(alpha-1)/(1-U_t)*u_t*(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)^2
    C_ss2=-alpha*(1-U_t_1)^(alpha-1)*u_t_1*(Y[1:n-1]-mu)*(Y[2:n]-mu)
    C_ss3=(1-(1-U_t_1)^alpha)*(Y[2:n]-mu)^3/sigma
    A_ss2=alpha*(1-U_t)^(alpha-1)*(C_ss1+C_ss2+C_ss3)*u_t/sigma^4
    A_ss=A_ss1+A_ss2
    
  }
  
  d2L_dmudsigma=function(mu,sigma,alpha){
    h=0.0000001
    (dL_dsigma(mu+h,sigma,alpha)-dL_dsigma(mu,sigma,alpha))/h
  }
  
  #d2L_dmudsigma(est[1],est[2],est[3])
  #-res$hessian[1,2]
  #h=0.0000001
    
  ############d^2L/dmudalpha#################################
  d2L_dmudalpha=function(mu,sigma,alpha){
    U_t_1=G(Y[1:n-1],mu,sigma);U_t=G(Y[2:n],mu,sigma)
    u_t_1=g(Y[1:n-1],mu,sigma);u_t=g(Y[2:n],mu,sigma)
    A=(1-U_t_1)^alpha+(1-U_t)^alpha-((1-U_t_1)^alpha)*((1-U_t)^alpha)
    A_a1=(1-U_t_1)^alpha*log(1-U_t_1)+(1-U_t)^alpha*log(1-U_t)
    A_a2=(1-U_t_1)^alpha*(1-U_t)^alpha*log((1-U_t_1)*(1-U_t))
    A_a=A_a1-A_a2
    A_mu1=alpha*(1-U_t_1)^(alpha-1)*(1-(1-U_t)^alpha)*u_t_1/sigma
    A_mu2=alpha*(1-U_t)^(alpha-1)*(1-(1-U_t_1)^alpha)*u_t/sigma
    A_mu=A_mu1+A_mu2
    
    B1=1-(1-U_t)^alpha+alpha*log(1-U_t_1)*(1-(1-U_t)^alpha)
    B2=-alpha*(1-U_t)^alpha*log(1-U_t)
    A_am1=(1-U_t_1)^(alpha-1)*(B1+B2)*u_t_1/sigma
    C1=1-(1-U_t_1)^alpha+alpha*log(1-U_t)*(1-(1-U_t_1)^alpha)
    C2=-alpha*(1-U_t_1)^alpha*log(1-U_t_1)
    A_am2=(1-U_t)^(alpha-1)*(C1+C2)*u_t/sigma
    A_am=A_am1+A_am2
    L1=A_am/(alpha-1+A)-(1+A_a)*A_mu/(alpha-1+A)^2
    L2=1/sigma*( u_t_1/(1-U_t_1)+u_t/(1-U_t) )-A_mu/alpha^2/A
    L3=(1/alpha-2)*( A_am/A-A_mu*A_a/A^2 )
    return((sum(L1)+sum(L2)+sum(L3))/n)
  }
  
  ############d^2L/dsigmadalpha (in progress) #################################
  d2L_dsigmadalpha=function(mu,sigma,alpha){
    h=0.0000001
    (dL_dsigma(mu,sigma,alpha+h)-dL_dsigma(mu,sigma,alpha))/h
  }
    
  #############Jacobian function######################################
  Ja=function(mu,sigma,alpha){
    AA=matrix( c(d2L_dmu2(mu,sigma,alpha),d2L_dmudsigma(mu,sigma,alpha),d2L_dmudalpha(mu,sigma,alpha),d2L_dmudsigma(mu,sigma,alpha),d2L_dsigma2(mu,sigma,alpha),d2L_dsigmadalpha(mu,sigma,alpha),d2L_dmudalpha(mu,sigma,alpha),d2L_dsigmadalpha(mu,sigma,alpha),d2L_dalpha2(mu,sigma,alpha)),3,3 )
    return(AA)
  }
  
  ###############Multivariate Newton Raphson#####################
  X=matrix(,1,3)
  tau=cor(Y[1:n-1],Y[2:n],method="kendall")
  alpha_est=-2*tau/(tau-1)
  alpha_est=2
  X[1,]=c(mean(Y),sd(Y),alpha_est)   #initial value
  i=2
  Ran.num=1
  
  repeat{
    Z=X
    X=matrix(,i,3)
    X[1:i-1,]=Z[1:i-1,]
    ##
    Aa=Ja(X[i-1,1],X[i-1,2],X[i-1,3])
    Ainv11=Aa[2,2]*Aa[3,3]-Aa[3,2]*Aa[2,3]
    Ainv12=Aa[1,2]*Aa[3,3]-Aa[3,2]*Aa[1,3]
    Ainv13=Aa[1,2]*Aa[2,3]-Aa[2,2]*Aa[1,3]
    Ainv21=Aa[2,1]*Aa[3,3]-Aa[3,1]*Aa[2,3]
    Ainv22=Aa[1,1]*Aa[3,3]-Aa[3,1]*Aa[1,3]
    Ainv23=Aa[1,1]*Aa[2,3]-Aa[1,3]*Aa[2,1]
    Ainv31=Aa[2,1]*Aa[3,2]-Aa[3,1]*Aa[2,2]
    Ainv32=Aa[1,1]*Aa[3,2]-Aa[1,2]*Aa[3,1]
    Ainv33=Aa[1,1]*Aa[2,2]-Aa[1,2]*Aa[2,1]
    Ainv=matrix(c(Ainv11,-Ainv21,Ainv31,-Ainv12,Ainv22,-Ainv32,Ainv13,-Ainv23,Ainv33),3,3)/det(Aa)
    ##
    X[i,]=X[i-1,]-Ainv%*%F(X[i-1,1],X[i-1,2],X[i-1,3])
    if(1*is.nan(X)[i,1]==1){
      X=matrix(,2,3)
      X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
      Ran.num=Ran.num+1
      i=1
    }else if(abs(X[i,1]-X[i-1,1])<0.0001&abs(X[i,2]-X[i-1,2])<0.0001&abs(X[i,3]-X[i-1,3])<0.0001&abs(X[i,3]-alpha_est)>5){
      X=matrix(,2,3)
      X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
      Ran.num=Ran.num+1
      i=1
    }else if(abs(X[i,1]-X[i-1,1])<0.0001&abs(X[i,2]-X[i-1,2])<0.0001&abs(X[i,3]-X[i-1,3])<0.0001&X[i,2]>0&abs(X[i,3]-alpha_est)<5){break
    }else if(abs(X[i,1]-X[i-1,1])<0.0001&abs(X[i,2]-X[i-1,2])<0.0001&abs(X[i,3]-X[i-1,3])<0.0001&X[i,2]<0){
      X=matrix(,2,3)
      X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
      Ran.num=Ran.num+1
      i=1
    }else if(abs(X[i,1]-X[i-1,1])>10^10&abs(X[i,2]-X[i-1,2])>10^10&abs(X[i,3]-X[i-1,3])>10^10){
      X=matrix(,2,3)
      X[1,]=c(mean(Y),sd(Y),alpha_est+runif(1,-D,D))   #initial value
      Ran.num=Ran.num+1
      i=1
    }
    if(Ran.num>=100){break}
    i=i+1
  }
  mle.res=X[length(X[,1]),]
  
  if(Ran.num>=10){ mle.res=c(mean(Y),sd(Y),alpha_est) }
  
  UCL=mle.res[1]+k*mle.res[2]
  LCL=mle.res[1]-k*mle.res[2]
  
  result=c(mu=mle.res[1],sigma=mle.res[2],alpha=mle.res[3],UCL=UCL,LCL=LCL)
  
  ####### Plot Control Chart #######
  if(plot==TRUE){
    Min=min(min(Y),LCL)
    Max=max(max(Y),UCL)
    ts.plot(Y,type="b",ylab="Y",ylim=c(Min,Max))
    abline(h=result[1])
    abline(h=UCL,lty="dotted",lwd=2)
    abline(h=LCL,lty="dotted",lwd=2)  
    text(0,LCL+(result[1]-LCL)*0.1,"LCL")
    text(0,UCL-(UCL-result[1])*0.1,"UCL")
  }
  
  out_control=which(  (Y<LCL)|(UCL<Y)  )
  if(length(out_control)==0){out_control="NONE"}
  
  Gradient=F(mle.res[1],mle.res[2],mle.res[3])
  Hessian=Ja(mle.res[1],mle.res[2],mle.res[3])
  
 list(estimates=result,out_of_control=out_control,
         Gradient=Gradient,Hessian=Hessian,Mineigenvalue_Hessian=min(eigen(Hessian)$value))
  
}
