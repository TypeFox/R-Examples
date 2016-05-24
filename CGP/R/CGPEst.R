CGPEst <-
function(X,yobs,nugget_l=0.001,num_starts=5,theta_l=NULL,alpha_l=NULL,kappa_u=NULL){
  
  DD<-as.matrix(X)
  var_names<-colnames(DD)
  
  n<-nrow(DD)
  p<-ncol(DD)
  one<-rep(1,n)
  onem<-rep(1,n-1)
  
  #Standardized design matrix
  Stand_DD<-apply(DD,2,function(x) (x-min(x))/max(x-min(x)) )
  
  #Scales for standardization
  scales<-apply(DD,2,function(x) max(x)-min(x))
  
  
  if(length(theta_l)>1) print("Error: Lower bound for theta needs to be specified as a scalar!")
  if(length(alpha_l)>1) print("Error: Lower bound for alpha needs to be specified as a scalar!")
  if(length(kappa_u)>1) print("Error: Upper bound for kappa needs to be specified as a scalar!")
  
  
  if(is.null(theta_l)) theta_l<-0.0001
  theta_lower<-rep(theta_l,p)
  
  #d_avg<-sqrt(1/mean(1/dist(Stand_DD)^2))
  if(is.null(alpha_l)) alpha_l<-log(10^2)*mean(1/dist(Stand_DD)^2)
  
  theta_upper<-rep(alpha_l,p)
  kappa_l<-alpha_l
  
  if(is.null(kappa_u)) kappa_u<-log(10^6)*mean(1/dist(Stand_DD)^2)
  
  if(sum(theta_l>alpha_l)>0) print("Error: Lower bound of theta exceeds the upper bound!")
  
  lower=c(nugget_l,theta_lower,kappa_l,0)
  upper=c(1,theta_upper,kappa_u,1)
  
  
  # Function to construct the (n by n) correlation matrix.
  PSI<-function(theta){
    A<-DD%*%diag(sqrt(theta),ncol=p)
    A<-as.matrix(dist(A, diag=T, upper=T))
    R<-exp(-A^2)
    return(R)
  }
  Stand_PSI<-function(Stand_theta){
    A<-Stand_DD%*%diag(sqrt(Stand_theta),ncol=p)
    A<-as.matrix(dist(A, diag=T, upper=T))
    R<-exp(-A^2)
    return(R)
  }
  
  # Likelihood function
  var.MLE.DK<-function(ww){
    lambda<-ww[1]
    Stand_theta<-ww[2:(p+1)]
    kappa<-ww[p+2]
    Stand_alpha<-kappa+Stand_theta
    bw<-ww[p+3]
    G<-Stand_PSI(Stand_theta)
    L<-Stand_PSI(Stand_alpha)
    Gbw<-Stand_PSI(Stand_theta*bw)
    Sig<-diag(n)
    for(rep in 1:4){
      Q<-G+lambda*Sig^(1/2)%*%L%*%Sig^(1/2)
      invQ <- solve(Q)
      beta <- (one %*% invQ %*% yobs)/(one %*% invQ %*% one)
      temp<-invQ%*%(yobs-beta*one)
      gip<-beta*one+G%*%temp
      e<-yobs-gip
      Sig<-diag(c(Gbw%*%e^2/(Gbw%*%one)))
      Sig2<-mean(diag(Sig))
      Sig<-Sig/Sig2
    }
    Q<-G+lambda*Sig^(1/2)%*%L%*%Sig^(1/2)
    invQ <- solve(Q)
    beta <- (one %*% invQ %*% yobs)/(one %*% invQ %*% one)
    tau2<- t(yobs-beta*one)%*%invQ%*%(yobs-beta*one)/n
    val<- log(det(Q))+ n*log(tau2) 
    if(!is.finite(val)) val<-1000000
    return(val)
  }
  
    
  # Optimize the likelihood function using "num_starts" Latin-hypercube sampled points as random starts. 
  n_par<-p+3
  n_candidate<-500+num_starts
  LHD<-function(N,k){
    x<-matrix(rep(1:N,k),ncol=k,nrow=N)
    x<-apply(x,2,sample)
    x<-(x-0.5)/N
    return(x)
  }
  starts<-LHD(n_candidate,n_par)
  #require(lhs)
  #starts<-maximinLHS(n_candidate,n_par)
  range_S<-matrix(rep(upper-lower,n_candidate),ncol=n_par,byrow=TRUE)
  low_S<-matrix(rep(lower,n_candidate),ncol=n_par,byrow=TRUE)
  starts<-starts*range_S+low_S
  cand_obj<-apply(starts,1,var.MLE.DK)
  index<- (rank(cand_obj,ties.method="min")<=num_starts)
  Starts<-starts[index,]
  Fop<-function(x) optim(x,var.MLE.DK,lower=lower,upper=upper,method="L-BFGS-B")$val
  op_obj<-apply(Starts,1,Fop)
  beststart<-Starts[which.min(op_obj),]
  op<-optim(beststart,var.MLE.DK,lower=lower,upper=upper,method="L-BFGS-B")
  
  objval<-op$val
  lambda<-op$par[1]
  Stand_theta<-op$par[2:(p+1)]
  kappa<-op$par[p+2]
  Stand_alpha<-kappa+Stand_theta
  bw<-op$par[p+3]
  
  #Adjust the scales of theta and alpha for unstandaridized designs DD. 
  theta<-Stand_theta/scales^2
  alpha<-Stand_alpha/scales^2
  
  
  ###########Leave-one-out Cross Validation Start###########################
  Yjfp<-rep(0,n)
  
  for(jf in 1:n){
    
    G<-PSI(theta)[-jf,-jf]
    L<-PSI(alpha)[-jf,-jf]
    Gbw<-PSI(theta*bw)[-jf,-jf]
    Sig<-diag(n-1)
    
    for(rep in 1:4){
      Q<-G+lambda*Sig^(1/2)%*%L%*%Sig^(1/2)
      invQ <- solve(Q)
      beta <- (onem %*% invQ %*% yobs[-jf])/(onem %*% invQ %*% onem)
      temp<-invQ%*%(yobs[-jf]-beta*onem)
      gip<-beta*onem+G%*%temp
      e<-yobs[-jf]-gip
      Sig<-diag(c(Gbw%*%e^2/(Gbw%*%onem)))
      Sig2<-mean(diag(Sig))
      Sig<-Sig/Sig2
      #print(diag(Sig))
    }
    
    Q<-G+lambda*Sig^(1/2)%*%L%*%Sig^(1/2)
    invQ <- solve(Q)
    beta <- (onem %*% invQ %*% yobs[-jf])/(onem %*% invQ %*% onem)
    tau2<- t(yobs[-jf]-beta*onem)%*%invQ%*%(yobs[-jf]-beta*onem)/(n-1)
    temp<-invQ%*%(yobs[-jf]-beta*onem) 
    
    g<-PSI(theta)[jf,-jf]
    l<-PSI(alpha)[jf,-jf]
    gbw<-PSI(theta*bw)[jf,-jf]
    
    vjf<-(t(gbw)%*%(e^2)/(t(gbw)%*%onem))/Sig2
    vjf<-as.vector(vjf)
    q<-g+lambda*sqrt(vjf)*Sig^(1/2)%*%l
    Yjfp[jf]<-beta+t(q)%*%temp
    
  }
  
  rmscv<- sqrt(sum((yobs-Yjfp)^2)/n)
  
  ###############Leave-one-out Cross Validation End##################
  
  
  G<-PSI(theta)
  L<-PSI(alpha)
  Gbw<-PSI(theta*bw)
  
  # Obtain the local volatility matrix "Sig"
  Sig<-diag(n) 
  for(rep in 1:4){
    Q<-G+lambda*Sig^(1/2)%*%L%*%Sig^(1/2)
    invQ <- solve(Q)
    beta <- (one %*% invQ %*% yobs)/(one %*% invQ %*% one)
    temp<-invQ%*%(yobs-beta*one)
    gip<-beta*one+G%*%temp
    e<-yobs-gip
    res2<-e^2
    Sig<-diag(c(Gbw%*%res2/(Gbw%*%one)))
    sf<-mean(diag(Sig))
    Sig<-Sig/sf
    #print(diag(Sig))
  }
  
  
  Q<-G+lambda*Sig^(1/2)%*%L%*%Sig^(1/2)
  invQ <- solve(Q)
  beta <- (one %*% invQ %*% yobs)/(one %*% invQ %*% one)
  tau2<- t(yobs-beta*one)%*%invQ%*%(yobs-beta*one)/n
  temp<-invQ%*%(yobs-beta*one)
  
  theta<-matrix(theta,nrow=1)
  alpha<-matrix(alpha,nrow=1)
  if(!is.null(var_names)){
    colnames(theta)=var_names
    colnames(alpha)=var_names
  }
  
  val<-list(X=DD,yobs=yobs,var_names=var_names,lambda=lambda,theta=theta,alpha=alpha,bandwidth=bw,Sig_matrix=Sig,sf=sf,res2=res2,temp_matrix=temp,invQ=invQ,mu=beta,tau2=tau2,beststart=beststart,objval=objval,rmscv=rmscv,Yp_jackknife=Yjfp)
  
  return(val)
}
