######################
# R libraries        #
######################
library(quadprog); library(quantreg)

#######################
#load R functions     # 
#######################

mspline_knots<-function(knots_inter,spline_df){
  n_knots<-length(knots_inter)+spline_df #number of basis functions
  knots_all<-c(rep(0,spline_df),knots_inter,rep(1,spline_df+1))
  return(knots_all)
}

ispline_knots<-function(knots_inter,spline_df){
  knots_all<-c(rep(0,spline_df),knots_inter,rep(1,spline_df+2))
  return(knots_all)
}

mspline<-function(tau, spline_df, m, M_knots){
  if ( tau<M_knots[m]  ||  tau>=M_knots[m+spline_df] ){v<-0} 
  else if (spline_df==1){v<- 1/(M_knots[m+1]-M_knots[m])}
  else{
    d1 <- tau-M_knots[m]
    d2 <- M_knots[m+spline_df]-tau
    v <- d1*mspline(tau, spline_df-1, m, M_knots) + d2*mspline(tau, spline_df-1, m+1, M_knots)
    v <- v*spline_df
    v <- v/((spline_df-1)*(M_knots[m+spline_df]-M_knots[m]))
  }
  return(v)
}

ispline<-function(tau, spline_df=3, m, I_knots){
  bin<-findInterval(tau,I_knots)  
  if (bin<m){v<-0}
  else if ( (bin-spline_df+1)>m ){v<-1}
  else{
    v <- 0
    for (l in m:(bin+1)){
      v <- v+mspline(tau, spline_df+1, l, I_knots)*(I_knots[l+spline_df+1]-I_knots[l])/(spline_df+1)
    }
  }
  return (v)
}

ispline2<-function(tau, spline_df, m, I_knots){
  bin<-findInterval(tau,I_knots)
  if(bin < m){v<-0}
  else if ((bin-spline_df+1)>m){v<-1}
  else if (bin == (m)){v<-   
    (tau - I_knots[m])^3/(
      (I_knots[m+1]  - I_knots[m])*(I_knots[m+2]  - I_knots[m])*(I_knots[m+3]  - I_knots[m])
      )
  }
  else if(bin == (m+1)){
    I1<-  (tau-I_knots[m])/(I_knots[m+2]-I_knots[m+1])
    I2<-  ((tau-I_knots[m+1])^2 - (I_knots[m+3]-tau)^2)/((I_knots[m+3]-I_knots[m+1])*(I_knots[m+2]-I_knots[m+1]))
    I3<-  (I_knots[m+3]-tau)^3/((I_knots[m+3]-I_knots[m+1])*(I_knots[m+3]-I_knots[m])*(I_knots[m+2]-I_knots[m+1]))-(tau-I_knots[m])^3/((I_knots[m+3]-I_knots[m])*(I_knots[m+2]-I_knots[m])*(I_knots[m+2]-I_knots[m+1]))
    v<- I1 + I2 + I3
      }
  else{v<- 
        1 - (I_knots[m+3]-tau)^3/((I_knots[m+3]-I_knots[m+2])*(I_knots[m+3]-I_knots[m+1])*(I_knots[m+3]-I_knots[m]))
        }  
    return(v)
}

mspline2<-function(tau, spline_df, m, I_knots){
  bin<-findInterval(tau,I_knots)
  if(bin < m){v<-0}
  else if ((bin-spline_df+1)>m){v<-0}
  else if (bin == (m)){v<-   
    3*(tau - I_knots[m])^2/((I_knots[m+1]-I_knots[m])*(I_knots[m+2]-I_knots[m])*(I_knots[m+3]-I_knots[m]))
  }
  else if(bin == (m+1)){
    I1<-  1/(I_knots[m+2]-I_knots[m+1])
    I2<-  (2*(I_knots[m+3]-I_knots[m+1]))/((I_knots[m+3]-I_knots[m+1])*(I_knots[m+2]-I_knots[m+1]))
    I3<-  -3 *(
      (I_knots[m+3]-tau)^2/((I_knots[m+3]-I_knots[m+1])*(I_knots[m+3]-I_knots[m])*(I_knots[m+2]-I_knots[m+1])) + (tau-I_knots[m])^2/((I_knots[m+3]-I_knots[m])*(I_knots[m+2]-I_knots[m])*(I_knots[m+2]-I_knots[m+1]))      
      )
    v<- I1 + I2 + I3
  }
  else{v<- 3*(I_knots[m+3]-tau)^2/((I_knots[m+3]-I_knots[m+2])*(I_knots[m+3]-I_knots[m+1])*(I_knots[m+3]-I_knots[m]))
  }  
  return(v)
}

make_M<-function(nu,tau,spline_df,M_knots){
  return(c(0,sapply(nu,mspline,spline_df = spline_df, M_knots = M_knots, tau=tau)))
}
make_I<-function(nu,tau,spline_df,I_knots){
  return(c(1,sapply(nu,ispline,spline_df = spline_df, I_knots = I_knots, tau=tau)))  
}

#function that creates my initial quantile estimates
quant_sim<-function(y,x,tau=seq(.1,.9,.1),crosscov=T){
  p<-2
  if(is.matrix(x)){p<-ncol(x)+1}#don't include an intercept in x!
  nt<-length(tau)
  #Do the individual regressions:
  beta<-matrix(0,p,nt)
  Hinv<-array(0,c(p,p,nt))
  for(t in 1:nt){
    fit<-rq(y~x,tau=tau[t])
    beta[,t]<-fit$coef
    Hinv[,,t]<-summary(fit,cov=T,se="ker")$Hinv
  }
  #put the pieces together:
  J_partial<-summary(rq(y~x,tau=0.5),cov=T,se="ker")$J
  beta_hat<-as.vector(t(beta))
  beta_cov<-matrix(0,p*nt,p*nt)
  taus<-kronecker(rep(1,p),tau)
  label<-kronecker(1:p-1,rep(1,nt))
  for(t1 in 1:nt){
    for(t2 in 1:nt){
      if(t1==t2 | crosscov){
        covar<-Hinv[,,t1]%*%J_partial%*%Hinv[,,t2]
        covar<-covar*(min(tau[c(t1,t2)])-tau[t1]*tau[t2])
        beta_cov[taus==tau[t1],taus==tau[t2]]<-covar
      }
    }
  }
  return(list(beta_hat=beta_hat,beta_cov=beta_cov,tau=taus,label=label))
}

quant_sim_tails<-function(y,x,tau,crosscov=T){
  p<-2
  if(is.matrix(x)){p<-ncol(x)+1}#don't include an intercept in x!
  beta<-matrix(0,p,1)
  Hinv<-array(0,c(p,p,1))
  fit<-rq(y~x,tau)
  beta_hat<-fit$coef
  beta_cov<-summary(fit,cov=T,se="iid")$cov
  return(list(beta_hat=beta_hat,beta_cov=beta_cov))
}

mle_start<-function(betahat,II,N_tau,P,M){
  
  #vector appearing in the quadratic function to be minimized (first element of normal equations)
  #W has N_tau*(1+G+P) rows
  
  W<- matrix(
    cbind(II,
          matrix(0,nrow=(nrow(II)),ncol=(2*ncol(II))))
    ,nrow=(nrow(II))
  )
  
  dvec<- matrix(t(betahat)%*%W,ncol=1)
  #second element of normal equations
  Dmat<- t(W)%*%W; diag(Dmat)[  (M*P+1):ncol(W)]<-  10^-5
  
  ##################    equality constraints ###################################
  
  #ensure the nonintercept alphas equal the difference of the betas and gammas
  #must hold for all m
  A1<- cbind(diag(P*M),-diag(P*M),diag(P*M))
  #################   inequality constraints ##################################
  #ensure the intercept is larger than the negative slopes
  #doesn't have to hold for m = 1
  A2<-kronecker(t(c(rep(0,P),1,rep(0,P-1),rep(-1,P))),diag(M))
  A2<-A2[2:M,]
  #pick up here#
  #ensure all betas and gammas are positive (which ensures all intercept alphas are positive)
  #must hold for all m
  A3<-diag(3*P*M)
  A3<-A3[(P*M+1):(3*P*M),]
  Amat<- matrix(rbind(A1,A2,A3),ncol=ncol(A1))
  bvec<- c(rep(0,nrow(A1)),      rep(.1, nrow(Amat)-nrow(A1)  ))
  #ncol Amat = length(dvec)
  theta_QP<-solve.QP(Dmat, dvec, t(Amat), bvec,meq=(nrow(A1)), factorized=FALSE)$solution
  return(theta_QP)
}
#new inequality constraints
mle_start_2<-function(betahat,II,N_tau,P,M){
  
  #vector appearing in the quadratic function to be minimized (first element of normal equations)
  #W has N_tau*(1+G+P) rows
  
  W<- matrix(
    cbind(II,
          matrix(0,nrow=(nrow(II)),ncol=(2*ncol(II))))
    ,nrow=(nrow(II))
    )
  
  dvec<- matrix(t(betahat)%*%W,ncol=1)
  #second element of normal equations
  Dmat<- t(W)%*%W; diag(Dmat)[  (M*P+1):ncol(W)]<-  10^-5
  
  ##################    equality constraints ###################################
  
  #ensure the alphas equal the difference of the betas and gammas
  #must hold for all m
  A1<- cbind(diag(P*M),-diag(P*M),diag(P*M))
  #################   inequality constraints ##################################
  #ensure the intercept is larger than the positive and the negative slopes
  #doesn't have to hold for m = 1
#   A2<-kronecker(t(c(1,rep(0,P-1),rep(-1,P),rep(-1,P))),diag(M))
#   A2<-A2[2:M,]
  A2<-kronecker(t(c(rep(0,P),1,rep(-1,P-1),rep(-1,P))),diag(M))
  A2<-A2[2:M,]
  #pick up here#
  #ensure all betas and gammas are positive (which ensures all intercept alphas are positive)
  #must hold for all m
  A3<-diag(3*P*M)
  A3<-A3[(P*M+1):(3*P*M),]
  Amat<- matrix(rbind(A1,A2,A3),ncol=ncol(A1))
  bvec<- c(rep(0,nrow(A1)),      rep(.01, nrow(Amat)-nrow(A1)  ))
  #ncol Amat = length(dvec)
  theta_QP<-solve.QP(Dmat, dvec, t(Amat), bvec,meq=(nrow(A1)), factorized=FALSE)$solution
  return(theta_QP)
}

#function that evaluates the quantile function for a given tau, x and alpha
#w is a vector of length M
Q<-function(tau,w, I_knots){
  M<-length(w)
  return(t(w)%*%make_I(nu = 1:(M-1),tau = tau, spline_df = 3,I_knots = I_knots))
}
#function that evaluates the quantile density for a given tau, x and alpha
#w is a vector of length (M-1)
q<-function(tau,w,M_knots){
  M<-length(w)
  return(t(w[2:M])%*%make_M(nu = 1:(M-1),tau = tau, spline_df = 3,M_knots = M_knots))
}


#function that ensures our tau is in the proper interval
rootfind<-function(M_knots_length, I_knots_length, M, df_spline, tau_scalar, w, y_scalar,  M_knots, I_knots, bin){
  f <-        Q(tau_scalar,w, I_knots);
  f_prime <-  q(tau_scalar,w[2:M], M_knots);
  diff <- abs(y_scalar-f);
  iter <- 0;
  while(diff>.0000001 && iter < 1000)  {
    iter<- iter+1;
    tau_scalar <- tau_scalar  - (f-y_scalar)/f_prime;
    if(tau_scalar<.05||tau_scalar>.95){tau_scalar<-runif(1,min=.05,max=.95)}
    #one find_interval update should be sufficient
    bin<-  findInterval(x=tau_scalar, vec = I_knots);
    f <- Q(tau_scalar,w, I_knots);
    f_prime <-  q(tau_scalar,w[2:M], M_knots);
    diff = abs(y_scalar-f);
  }
  #ensure final solution is proper
  if(0.05 > tau_scalar || tau_scalar > 0.95){
    tau_scalar = runif(1, min = 0.05, max = 0.95);
    rootfind (M_knots_length, I_knots_length, M, df_spline, tau_scalar, w, y_scalar,  M_knots, I_knots, bin);
  }
  list(tau=tau_scalar,bin=bin)
}

pareto<-function(sigma, xi, v){
  return (log(0.05) -log(sigma) -(1/xi+1)*log(1+xi*v/sigma))
}

pareto_tau_low<-function(sigma, xi, v){
  return (  .05- .05*( 1- (1+xi*v/sigma)^(-1/xi))
  )
}
pareto_tau_high<-function(sigma, xi, v){
  return (  .95+ .05*( 1- (1+xi*v/sigma)^(-1/xi)))
}

exponential<-function(sigma, v){
  return (log(.05)-log(sigma) -(v/sigma))
}

exponential_tau_low<-function(sigma, v){
  return (  .05- .05*(1- exp(-v/sigma)));
}

exponential_tau_high<-function(sigma, v){
  return (  .95+ .05*(1- exp(-v/sigma)))
}

#build the covariance matrix
#sigma is the standard deviation
#rho is the correlation - must be less than sigma

make_AR<-function(rho,G){
  times<-1:G
  H <- abs(outer(times, times, "-"))
  return(rho^H)
  }

#returns the lowertriangular portion of a G x G AR-1 precision matrix with correlation rho
AR_chol<- function(rho,G){
  rho_inv<- 1/rho
  G_chol<- matrix(0,G,G)
  for(g in 1:G){G_chol[g,1]<- 1/(rho_inv^(g-1))}
  for(g in 2:G){G_chol[g,2]<- sqrt(rho_inv^2 -1)*G_chol[g,1]}
  for(g in 3:G){
    for(l in 3:g){
      G_chol[g,l]<- rho_inv*G_chol[g,l-1]
    }
  }
  for(g in 2:G){G_chol[g,g]<- sqrt(1-rho^2)}  
  return(G_chol)
}
#prior function for updating theta 
#PREC is G x G x M x P dimensional array

inverse<-function(OMEGA){
  return(chol2inv(chol(OMEGA))) 
}

sumprec<-function(rho,G){
	return(2*1/(1-rho^2)+(G-2)*(1+2*rho^2/(1-rho^2))+2*(G-1)*(-rho/(1-rho^2)))
}

qlogistic<-function(tau){log(tau/(1-tau))}

make.B<-function(tau,kappa,base,shape){

  L<-length(kappa)-1

  if(base == "Gaussian"){
    qt<-qnorm(tau)
    qk<-qnorm(kappa)
  }
  if(base == "t"){
    qt<-qt(tau,df=shape)
    qk<-qt(kappa,df=shape)
  }
  if(base == "logistic"){
    qt<-qlogistic(tau)
    qk<-qlogistic(kappa)
  }
  if(base == "ALAP"){
    library(VGAM)
    qt<-qalap(tau,tau=shape)
    qk<-qalap(kappa,tau=shape)
  }
  if(base == "weibull"){
    qt<-qweibull(tau,shape=shape)
    qk<-qweibull(kappa,shape=shape)
  }
  if(base == "gamma"){
    qt<-qgamma(tau,shape)
    qk<-qgamma(kappa,shape)
  }

  B<-matrix(0,length(tau),L)
  for(l in 1:L){
    low<-kappa[l]
    high<-kappa[l+1]
    if(l==1){
      q<-ifelse(tau<=high,qt,qk[2])
    }
    if(l>1){
      q<-ifelse(tau<=low,0,qt-qk[l])
      q<-ifelse(tau>high,qk[l+1]-qk[l],q)
    }
    B[,l]<-q
  }
return(B)}

log_pdf<-function(y,X,beta,alpha_star,kappa,base,shape){

   L<-length(kappa)-1
   n<-length(y)

   if(base == "Gaussian"){
     q0<-function(tau,mn,s,shape){qnorm(tau,mn,s)}
     d0<-function(y,mn,s,shape){dnorm(y,mn,s,log=T)}
   }
   if(base == "t"){
     q0<-function(tau,mn,s,shape){mn+s*qt(tau,df=shape)}
     d0<-function(y,mn,s,shape){dt((y-mn)/s,df=shape,log=T)-log(s)}
   }
   if(base == "ALAP"){
     library(VGAM)
     q0<-function(tau,mn,s,shape){mn+s*qalap(tau,tau=shape)}
     d0<-function(y,mn,s,shape){dalap((y-mn)/s,tau=shape,log=T)-log(s)}
   }
   if(base == "logistic"){
     q0<-function(tau,mn,s,shape){mn+s*log(tau/(1-tau))}
     d0<-function(y,mn,s,shape){
        z<-(y-mn)/s
        return(-z-2*log(1+exp(-z))-log(s))}
   }
   if(base == "weibull"){
     q0<-function(tau,mn,s,shape){qweibull(tau,shape=shape,scale=s)+mn}
     d0<-function(y,mn,s,shape){dweibull(y-mn,shape=shape,scale=s,log=T)}
   }
   if(base == "gamma"){
     q0<-function(tau,mn,s,shape){qgamma(tau,shape,1/s)+mn}
     d0<-function(y,mn,s,shape){dgamma(y-mn,shape,1/s,log=T)}
   }

   # truncation to ensure a valid quantile process
   alpha<-alpha_star
   for(l in 1:L){
     alpha[,l]<-threshold(alpha_star[,l])
   }
 
   breaks<-t(q_split(kappa,kappa,X,beta,alpha,base=base,shape=shape))
   bin<-rep(1,n)
   for(l in 2:L){
      bin<-ifelse(y>breaks[,l],l,bin)
   }


   s<-X%*%alpha
   s<-s[cbind(1:n,bin)]
   
   m1<-X%*%beta
   m2<-breaks[cbind(1:n,bin)]
   m3<-s*q0(kappa[bin],0,1,shape=shape)
   mn<-ifelse(bin==1,m1,m2-m3)
 
   lll<-d0(y,mn=mn,s=s,shape=shape)

return(lll)}

q_split<-function(tau,kappa,X,beta,alpha,base,shape){

    B<-make.B(tau,kappa,base=base,shape=shape)
    loc<-X%*%beta
    scale<-X%*%alpha
    Q<-matrix(loc,length(tau),nrow(X),byrow=T)+
       B%*%t(scale)

return(Q)}

threshold<-function(alpha_star,epsilon=0.001){
    a<-alpha_star
    negparts<-ifelse(alpha_star[-1]>0,0,alpha_star[-1])
    worst_case<-alpha_star[1]+sum(negparts)
    if(worst_case<epsilon){
       a[1]<-a[1]-worst_case+epsilon
    }
return(a)}



