qreg_spline <-
function(X,Y=NULL,Y_low=NULL,Y_high=NULL,status = NULL,
                  knots_inter = c(.1,.5,.9),         
                  Pareto = TRUE,
                  varying_effect=NULL,
                  tau= seq(0.05,0.95,0.05),
                  burn=10000, iters=50000,
                  q_low = 0.01, q_high = 0.99, 
                  sig_a = .1, sig_b = .1,
                  mu_var = 10^2, cbf_var = 10^3, tail_mean = -1.2, tail_var = .4,
                  cbf_eps = 0.5, theta_eps = 0.5,
                  tuning_theta = 1, tuning_tail = rep(1,2),
                  cred_low = .025, cred_high = .975,
                  seed = 1, verbose = TRUE
                  ){
  library(quantreg); library(quadprog)
  n1 <- n<- length(Y)
  n2<-n3<-n4<-n5<-0
  
  
  if(length(status)>0){
    
    
    n1<- sum(status==0)
    n2<- sum(status==1)
    n3<- sum(status==2)
    n4<- sum(status==3)
    n5<- sum(status==4)
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
    print("Invalid number of covariates affecting the scale")
    stop()
  }
  
  if(sum(is.na(X))>0){
    print("Must remove observations with missing covariates")
    stop()
  }
  if(mean(X[,1]==1)<1){
    print("First column of X must be all ones")
    stop()
  }
  if(min(abs(X))>1){
    print("All values of X much be between -1 and 1")
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
  if(q_low > .1){
    print("Lower threshold cannot exceed 0.1")
    stop()
  }
  if(q_high < .9){
    print("Upper threshold must be at least 0.9")
    stop()
  }
  if(sig_a < 0 || sig_b < 0 || mu_var < 0 || cbf_var < 0 || tail_var < 0){
    print("All scale, shape and variance hyperparameters must be positive")
  }
  if(cbf_eps<.1 || cbf_eps > .5 || theta_eps < .1 || theta_eps > .5){
    print("All resolvant parameters must be in [0.1,0.5]")
    stop()
  }
  if(tuning_theta < 0 || min(tuning_tail) < 0){
    print("All candidate variances must be positive")
    stop()
  }
  burn<- as.integer(burn)
  sweeps<- as.integer(iters)
  
  if(burn < 1 || sweeps < 1){
    print("Need positive integer for burns and iters")
    stop
  }
  
  if(cred_low < 0){
    print("need lower limits of credible interval in (0,1)")
    stop
  }
  if(cred_high > 1){
    print("need upper limits of credible interval in (0,1)")
    stop
  }
  if(cred_low >= cred_high){
    print("need cred_low < cred_high")
    stop
  }
  
  if(length(knots_inter)> length(unique(knots_inter))){
    print("knots must be unique")
    stop
  }
  if(q_low <= 0){
    print("lower threshold must be greater than 0")
    stop
  }
  if(q_high >= 1){
    print("upper threshold must be less than 1")
    stop
  }
  
  

  
  xi_zero<- 1 - Pareto
  tau<-sort(tau)
  knots_inter<-sort(knots_inter)

  N_tau<-length(tau)
  
  tau_start <- seq(.05,.95,.05)
  N_tau_start<- length(tau_start)
  
  N<-nrow(X)  
  
  spline_df<-3
  M<- 1+ spline_df+length(knots_inter) # number of spline basis functions plus intercept

  
  
  M_knots<-mspline_knots(knots_inter,spline_df)
  I_knots<-ispline_knots(knots_inter,spline_df)
  
  M_knots_length<-length(M_knots)
  I_knots_length<-length(I_knots)
  
  III_start<-matrix(0,nrow=length(tau_start),ncol=length(knots_inter)+spline_df)
  for(i in 1:ncol(III_start)){
    III_start[,i]<-t(sapply(tau_start,ispline,spline_df=spline_df,m=i,I_knots=I_knots))  #constructs the spline basis evaluated at the Tau1 knots
  }
  III_start<-cbind(1,III_start)  
  II<-kronecker(diag(P),III_start)
  
  I_low<- make_I(nu=1:(M-1),tau = q_low, spline_df = 3,  I_knots = I_knots)
  I_high<- make_I(nu=1:(M-1),tau = q_high, spline_df = 3,  I_knots = I_knots)
  M_low<- make_M(nu=1:(M-1),tau=q_low,spline_df = 3, M_knots = M_knots) 
  M_high<- make_M(nu=1:(M-1),tau=q_high,spline_df = 3, M_knots = M_knots)
  M_low <- q_low * M_low
  M_high <-  (1 - q_high) * M_high

  dummy_y <- rep(0,N)
  if(length(status)==0){status<- rep(0,N)}
  if(sum(status==0) > 0){dummy_y[c(status==0)]<-Y[c(status==0)]}
  if(sum(status==1) > 0){set.seed(seed); dummy_y[c(status==1)]<- Y[c(status==1)] + rnorm(sum(status==1),.1)  } #left censored
  if(sum(status==2) > 0){set.seed(seed); dummy_y[c(status==2)]<- Y[c(status==2)] + rnorm(sum(status==2),.1)  } #right censored
  if(sum(status==3) > 0){set.seed(seed); dummy_y[c(status==3)]<- (Y_low[c(status==3)] + Y_high[c(status==3)])/2 + rnorm(sum(status==3),.1)  } #interval censored
  
  fit<-0

  betahat_matrix<-matrix(0,N_tau_start*P,1)
suppressWarnings(    tryCatch(
          fit<-quant_sim(dummy_y,X[,2:P],tau=tau_start),
             error=function(e){print("Initial fit failed");return(fit)}
      )
)
  if(length(fit)==1){
    fit<- list(beta_hat=rep(1,N_tau_start*P))
  }

  betahat_matrix[,1]<-fit$beta_hat #Koenker's quantile estimates
  
  theta_start<- array(0,dim=c(M,P,3))
  theta_start[1,1,1]<- min(dummy_y)
  theta_start[2:M,1,1]<- (max(dummy_y) - min(dummy_y))/(P-1)
  
  if(var(betahat_matrix) !=0){
  tryCatch(
    theta_start[,,]<-mle_start_2(c(betahat_matrix[,1]),II,N_tau_start,P,M),
    error=function(e){return(theta_start)})
  }
  thetastar<-theta_start[,,1]
  if(P1 < P){thetastar[,(P1 + 1):P] <- 0}

  #starting value for sampler
  old_tau<-rep(.5,N)
  bin<-rep(spline_df+2,N)
  iter_flag<-0
  tuning_parms = array(10,dim=c(M,P))
  
  if(!is.null(Y)){Y[is.na(Y)]<-0}
  if(!is.null(Y_low)){Y_low[is.na(Y_low)]<-0}
  if(!is.null(Y_high)){Y_high[is.na(Y_high)]<-0}

  M_1 <- M - 1
  IKM<- matrix(0,M_1,6)
  for(m in 1:(M-1)){
    IKM[m,1]  <- (I_knots[m+1]  - I_knots[m]) * (I_knots[m + 2]  - I_knots[m]) * (I_knots[m + 3]  - I_knots[m])
    IKM[m,2]  <- (I_knots[m+2] - I_knots[m+1])
    IKM[m,3]  <- ((I_knots[m+3] - I_knots[m+1]) * (I_knots[m+2] - I_knots[m+1])) 
    IKM[m,4]  <- ((I_knots[m+3] - I_knots[m+1]) * (I_knots[m+3] - I_knots[m]) * (I_knots[m+2]-I_knots[m+1]))
    IKM[m,5]  <- -((I_knots[m+3] - I_knots[m]) * (I_knots[m+2] - I_knots[m]) * (I_knots[m+2] - I_knots[m+1]))
    IKM[m,6]  <- -((I_knots[m+3] - I_knots[m+2]) * (I_knots[m+3] - I_knots[m+1]) * (I_knots[m+3] - I_knots[m]))
  }
  IKM<- ifelse(IKM==0,0,IKM^-1)
  
  MKM<- IKM
  MKM[,1]<- 3*MKM[,1]
  MKM[,2]<- 3*IKM[,2]
  MKM[,4]<- -3 * MKM[,4]
  MKM[,5]<- 3 * MKM[,5]
  MKM[,6]<- -3 * MKM[,6]
  
  thetastar<- array(0,dim=c(M,P))
  thetastar[,1]<-1
  
  theta_keep<- array(dim=c(iters, M * P, 1))
  tuning_parms_keep<- array(dim=c(M * P,1, 1))
  acc_theta_keep <- array(dim = c(M *  P, 1))
  
  ############## IP prior inputs #######################
  
  mu<- array(0,dim=c(1,P)) 
  sigma2 <-  array(1,dim=c(1,P))
  rho <- array(.5,dim=c(1,P))
  mu_keep <- array(dim=c(iters, 1 *  P, 1)) 
  sigma2_keep<-rho_keep<-array(dim=c(iters,1*P,1))  
  
  ############### tail inputs #################
  xi_low <- xi_high <- 1
  
  set.seed(seed)
  
  
  tick<- proc.time()
  out<- MCMC_IP_C(burn,iters,
                  tuning_parms, tuning_tail,
                  cbf_eps, theta_eps, 
                  M,  P,  P1, N,  n1, n2, n3, n4, n5,
                  Y,  Y_low, Y_high, X,
                  M_knots_length, I_knots_length, spline_df, 
                  M_knots, I_knots,
                  M_low,  M_high, I_low,  I_high,   
                  thetastar,  mu,  sigma2,  rho,  xi_low,  xi_high,
                  q_low,  q_high, xi_zero,
                  mu_var, cbf_var, tail_mean, tail_var,
                  sig_a, sig_b,
                  IKM, MKM, M_1,verbose)
  tock<- proc.time()
  MCMC_time <- tock - tick
  if(verbose){print("Compiling results...")}
  
  # common values
  post_theta<- matrix(out$THETA,nrow=out$iters,byrow=T)
  tuning_parms<- out$tuning_parms
  acc_theta <- out$ACC_THETA/out$ATT_THETA  
  theta_out<- round(out$THETA_OUT/out$sweeps,1) #n_times theta was outside of parameter space
  theta_in<-  round(out$THETA_IN/out$sweeps,1)
  
  # prior values
  post_mu<- matrix(out$MU,nrow=out$iters,byrow=T)
  post_sigma2<- matrix(out$SIGMA2,nrow=out$iters,byrow=T)
  post_rho<- matrix(out$RHO,nrow=out$iters,byrow=T)
  rho_accepted<-   round(matrix(out$ACC_RHO/out$iters,1,nrow=P),2)
  
  # tail values
  post_xi_low <- out$XI_LOW
  post_xi_high <- out$XI_HIGH
  
  ############################# Part 2: Credible Intervals
  post_q<- array(0,dim=c(iters,N_tau,ncol=P))
  post_q_lower<- post_q_upper <- post_q_mean<- matrix(0,N_tau,P)

  I_high<- make_I(nu=1:(M-1),tau = q_high, spline_df = 3,  I_knots = I_knots)
  M_low<- make_M(nu=1:(M-1),tau=q_low,spline_df = 3, M_knots = M_knots) 
  M_high<- make_M(nu=1:(M-1),tau=q_high,spline_df = 3, M_knots = M_knots)

  N_tau_low <- sum(tau < q_low)
  N_tau_mid <- sum( (q_low <=   tau)*(tau  <= q_high))
  N_tau_high <- sum(tau > q_high)

  if(N_tau_low > 0){
    if(Pareto){ 
        for(i in 1:N_tau_low){
        I_low<- make_I(nu=1:(M-1),tau = q_low, spline_df = 3,  I_knots = I_knots) 
        M_low<- make_M(nu=1:(M-1),tau=q_low,spline_df = 3, M_knots = M_knots) 
        for(p in 1:P){
          low_scale <- post_theta[,((p-1)*M+1) : (p*M)]%*%M_low
          Q_low <- post_theta[,((p-1)*M+1) : (p*M)]%*%I_low
          post_q[,i,p] <-  Q_low - (q_low/post_xi_low) * low_scale  * ((tau[i] / q_low)^-post_xi_low -1)        
        }
      }
    }
    else{ 
      for(i in 1:N_tau_low){
        I_low<- make_I(nu=1:(M-1),tau = q_low, spline_df = 3,  I_knots = I_knots) 
        M_low<- make_M(nu=1:(M-1),tau=q_low,spline_df = 3, M_knots = M_knots) 
        low_scale <- post_theta[,1:M]%*%M_low
        Q_low <- post_theta[,1:M]%*%I_low
        post_q[,i,1] <-  Q_low +  q_low * low_scale  * log(tau[i]/q_low)        
        for(p in 2:P){
          low_scale <- post_theta[,((p-1)*M+1) : (p*M)]%*%M_low
          Q_low <- post_theta[,((p-1)*M+1) : (p*M)]%*%I_low
          post_q[,i,p] <-  Q_low +  q_low * low_scale  * log(tau[i]/q_low)        
        }
      }
    }
  }
  if(N_tau_mid > 0){
    tau_mid <- tau[(N_tau_low + 1):(N_tau_low + N_tau_mid)]
    III<-matrix(0,nrow=N_tau_mid,ncol=length(knots_inter)+spline_df)
    for(i in 1:ncol(III)){
      III[,i]<-t(sapply(tau_mid,ispline,spline_df=spline_df,m=i,I_knots=I_knots))  #constructs the spline basis evaluated at the Tau1 knots
    }
    III<-cbind(1,III)  
    for(p in 1:P){
      post_q[,(N_tau_low + 1):(N_tau_low + N_tau_mid),p]<- post_theta[,((p-1)*M+1) : (p*M)]%*%t(III)
    }
  }
  if(N_tau_high > 0){
    if(Pareto){ 
      for(i in (N_tau_low + N_tau_mid + 1):N_tau){
        I_high<-  make_I(nu=1:(M-1),tau = q_high, spline_df = 3,  I_knots = I_knots) 
        M_high<- make_M(nu=1:(M-1),tau=q_high,spline_df = 3, M_knots = M_knots) 
        for(p in 1:P){
          high_scale <- post_theta[,((p-1)*M+1) : (p*M)]%*%M_high
          Q_high <- post_theta[,((p-1)*M+1) : (p*M)]%*%I_high
          post_q[,i,p] <-  Q_high + (1 - q_high) * high_scale / post_xi_high   * ((( 1 - tau[i]) /(1 - q_high))^-post_xi_high -1)     
        }
      }
    }
    else{ 
      for(i in (N_tau_low + N_tau_mid + 1):N_tau){
        I_high<- make_I(nu=1:(M-1),tau = q_high, spline_df = 3,  I_knots = I_knots) 
        M_high<- make_M(nu=1:(M-1),tau=q_high,spline_df = 3, M_knots = M_knots) 
        high_scale <- post_theta[,1:M]%*%M_high
        Q_high <- post_theta[,1:M]%*%I_high
        post_q[,i,1] <-  Q_high +  q_high * high_scale  * log((1 - tau[i])/(1 - q_high))        
        for(p in 1:P){
          high_scale <- post_theta[,((p-1)*M+1) : (p*M)]%*%M_high
          Q_high <- post_theta[,((p-1)*M+1) : (p*M)]%*%I_high
          post_q[,i,p] <-  Q_high -  q_high * high_scale  * log((1 - tau[i])/(1 - q_high))        
        }
      }
    }
  }
  for(p in 1:P){
    post_q_lower[,p] <- apply(post_q[,,p],2,quantile,probs=cred_low)
    post_q_upper[,p]<- apply(post_q[,,p],2,quantile,probs=cred_high)
    post_q_mean[,p]<- apply(post_q[,,p],2,mean)
  }
  list(q=post_q,
       tau=tau,
       theta=post_theta,
       tuning_parms=tuning_parms, acc_theta = acc_theta, N_tau = N_tau,
       post_mu=post_mu, post_sigma2=post_sigma2, rho_keep=rho_keep, 
       post_xi_low = post_xi_low, post_xi_high = post_xi_high,
       q = post_q, q_lower = post_q_lower, q_upper = post_q_upper, q_mean = post_q_mean, 
       MCMC_time = MCMC_time,
       tau = tau,
       LPML=out$LPML,CPO=out$CPO,
       MCMC_time=MCMC_time,
       LPML = out$LPML,CPO=out$CPO,
       iters=iters,burn=burn
       )}
