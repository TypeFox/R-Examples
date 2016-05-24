qreg <-
function(X,Y=NULL,Y_low=NULL,Y_high=NULL,status = NULL,
         L=4,base="Gaussian",
         varying_effect=NULL,
         tau=seq(0.05,0.95,0.05),
         burn=10000, iters=50000){


beta_var = 100
shape_var = 1
mu_var = 100 
sig_a = 1
sig_b = 1
beta_can = 0.1
alpha_can = 0.1
shape_can = .1
beta_eps = 0.5
alpha_eps = 0.5
verbose=FALSE


  n1 <- n<- length(Y)
  n2<-n3<-n4<-0


  if(length(status)>0){

    
    n1<- sum(status==0)
    n2<- sum(status==1)
    n3<- sum(status==2)
    n4<- sum(status==3)
    n<- n1+n2+n3+n4    
    
    if(n != length(status)){
       print("Invalid entry for status")
       stop()
    }

    order<-c(which(status==0),which(status==1),which(status==2),which(status==3))
    if(!is.null(Y_low)){ Y_low<- Y_low[order]}
    if(!is.null(Y_high)){Y_high<- Y_high[order]}
    if(!is.null(Y)){     Y<-Y[order]}
    X<-X[order,]
  }

  N<-nrow(X)
  P<-ncol(X)
  P1<-varying_effect
  if(is.null(P1)){P1<-P}
  if(P1<1 || P1>P){
    print("Invalid number of covaraites affecting the scale")
    stop()
  }

  if(sum(is.na(X))>0){
    print("Must remove observations with missing covariates!")
    stop()
  }
  if(mean(X[,1]==1)<1){
    print("First column of X must be all ones!")
    stop()
  }
  if(min(abs(X))>1){
    print("All values of X much be between -1 and 1!")
    stop()
  }

  Yinit<-rowMeans(cbind(Y,Y_low,Y_high),na.rm=TRUE)

  sumY<-sum(Yinit)
  if(is.na(sumY)){
    print("Can't have missing responses")
  }
  if(abs(sumY)==Inf){
    print("Responses must be finite")
  }

  init_alpha<-matrix(0,P,L)
  if(base=="gamma" | base=="weibull"){
    init_beta<-rep(0,P)
    init_alpha[1,]<-sd(Yinit)
  }
  if(base!="gamma" & base!="weibull"){
    inits<-lm(Yinit~X-1)
    init_beta<-inits$coef
    init_alpha[1,]<-1.5*sd(inits$residuals)
    init_shape<-0.5
  }
  beta_can = beta_can*sd(Yinit)
  alpha_can = alpha_can*sd(Yinit)

  if(!is.null(Y)){Y[is.na(Y)]<-0}
  if(!is.null(Y_low)){Y_low[is.na(Y_low)]<-0}
  if(!is.null(Y_high)){Y_high[is.na(Y_high)]<-0}
  
  tuning_alpha <-rep(alpha_can,length(init_alpha))
  tuning_beta <- rep(beta_can,length(init_beta))
  tuning_shape <- shape_can

  tick<- proc.time()
  out<- MCMC_C(burn=burn, sweeps=iters, 
           beta=init_beta, alphastar=t(init_alpha), shape=0.5,
           X=X, y=Y, y_low=Y_low, y_high=Y_high, 
           tuning_alpha, tuning_beta, tuning_shape, 
           beta_eps=beta_eps, alpha_eps=alpha_eps,
           base=base, L=L, P=P, P1=P1, 
           N = N, n1 = n1, n2 = n2, n3 = n3, n4 = n4, 
           mu = rep(0, P), sigma2 = rep(1,P), rho = rep(.5, P), 
           beta_var=beta_var, shape_var=shape_var, mu_var=mu_var, 
           sig_a=sig_a, sig_b=sig_b, 
           verbose)
  tock<- proc.time()
  MCMC_time <- tock - tick
  if(verbose){print("Compiling results...")}

  acc_rate_alpha<-out$ACC_ALPHA/out$ATT_ALPHA
  acc_rate_beta<-out$ACC_BETA/out$ATT_BETA

  post_beta <-matrix(out$BETA,nrow=out$sweeps,byrow=T)
  post_alpha <-matrix(out$ALPHA,nrow=out$sweeps,byrow=T)
  post_shape <- matrix(out$SHAPE,nrow=out$sweeps,byrow=T)
  kappa<-seq(0,1,length=L+1)

  post_q<-array(0,c(out$sweeps,length(tau),P))
  for(j in 1:P){for(l in 1:length(tau)){
    post_q[,l,j]<-post_beta[,j]
  }}


  if(var(post_shape)==0){
     B<-make.B(tau,kappa,base,shape=out$post_shape[1])  
     for(j in 1:P1){
        these<-1:L + (j-1)*L
        post_q[,,j]<-post_q[,,j]+(post_alpha[,these])%*%t(B)
     }
  }

  if(var(post_shape)>0){
    for(i in 1:out$sweeps){ 
     B<-make.B(tau,kappa,base,shape=post_shape[i])  
     for(j in 1:P1){
        these<-1:L + (j-1)*L
        post_q[i,,j]<-post_q[i,,j]+B%*%(post_alpha[i,these])
     }
    }
  }

list(q=post_q,
     tau=tau,
     shape=post_shape,
     beta=post_beta,
     alpha=post_alpha,
     MCMC_time=MCMC_time,
     LPML = out$LPML,
     L=L,base=base,  kappa=kappa,
     iters=iters,burn=burn,
     acc_rate_alpha=out$ACC_ALPHA/out$ATT_ALPHA,
     acc_rate_beta=out$ACC_BETA/out$ATT_BETA
)}
