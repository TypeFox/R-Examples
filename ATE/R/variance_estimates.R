###############################################################################
#We begin with the estimating equations for the simple estimate
get.gk.simple<- function(X.vec,Y.scaler,Ti.int,lam,beta,tau1,tau0,
                         FUNu, rho1, N, n1, ...){
  lam<- as.vector(lam)
  be<- as.vector(beta)
  uk <- FUNu(as.numeric(X.vec))
  
  gk1<- Ti.int*rho1(crossprod(lam, uk),...)*uk - uk
  gk2<- (1-Ti.int)*rho1(crossprod(be, uk),...)*uk - uk
  gk3<- Ti.int*rho1(crossprod(lam,uk),...)*Y.scaler- tau1
  gk4<- (1-Ti.int)*rho1(crossprod(be, uk),...)*Y.scaler - tau0
  
  gk<- c(gk1,gk2,gk3,gk4) 
  gk
}

###############################################################################
#Secondly we have the estimating equation for the case of treatement of treated
get.gk.ATT<- function(X.vec,Y.scaler,Ti.int,beta,tau1,tau0,
                      FUNu, rho1,N,n1, ...){
  be<- as.vector(beta)
  uk <-  FUNu(as.numeric(X.vec))
  frac<-n1/N
  
  gk1<- (1-Ti.int)*rho1(crossprod(be, uk),...)*uk - 1/frac*Ti.int*uk
  gk2<- N/n1*(Ti.int*(Y.scaler- tau1))
  gk3<- (1-Ti.int)*rho1(crossprod(be, uk),...)*Y.scaler - tau0
  gk4<- Ti.int-frac
  
  gk<- c(gk1,gk2,gk3,gk4) 
  gk
  
}

###############################################################################
#Thirdly we have the estimating equation for the case of multiple treatment effects
#In this case we do not have a tau quantity, if we wish to compare two different 
#groups we can then formulate a variance estimate
get.gk.MT<- function(X.vec,Y.scaler,Ti.int,lam,
                     FUNu, rho1, big.J, taus, ...){
  lam<- as.vector(lam)
  uk <-  FUNu(as.numeric(X.vec))
  K<- length(lam)
  gk1<- rep(-uk,big.J)
  gk1[(Ti.int*K + 1):((Ti.int+1)*K)]<- rho1(crossprod(lam, uk),...)*uk + gk1[(Ti.int*K + 1):((Ti.int+1)*K)]
  gk2<- -taus
  gk2[Ti.int+1]<- gk2[Ti.int+ 1] + rho1(crossprod(lam, uk),...)*Y.scaler
  c(gk1,gk2)
}

###############################################################################
#This function obtains the big covariance matrix for all coefficients
#THe bottom right element of this will be the variance estimate for tau-simple
get.cov.simple<- function(X, Y, Ti, FUNu, rho, rho1, rho2, obj,...){
  N<- length(Y)
  n1<- sum(Ti)
  lam<- as.vector(obj$lam.p)
  be<- as.vector(obj$lam.q)
  tau1<- obj$Y1
  tau0<- obj$Y0
  umat <- t(apply(X, 1, FUNu))
  K<- length(lam)
  
  A<- matrix(0,ncol = 2*K, nrow = 2*K)
  C<- matrix(0, ncol = 2*K, nrow = 2)
  meat<- matrix(0, ncol = 2*K+2, nrow = 2*K+2)
  for(i in 1:N){
    A[(1:K),(1:K)]<- A[(1:K),(1:K)] + 
      Ti[i]*rho2(crossprod(lam, umat[i,]), ...)*tcrossprod(umat[i,])
    A[((K+1):(2*K)),((K+1):(2*K))]<- A[((K+1):(2*K)),((K+1):(2*K))] + 
      (1-Ti[i])*rho2(crossprod(be, umat[i,]),...)*tcrossprod(umat[i,])
    
    C[1,(1:K)]<- C[1,(1:K)] + Ti[i]*rho2(crossprod(lam, umat[i,]), ...)*Y[i]*umat[i,]
    C[2,((K+1):(2*K))]<- C[2,((K+1):(2*K))] +
      (1-Ti[i])*rho2(crossprod(be, umat[i,]), ...)*Y[i]*umat[i,]
    
    meat<- meat +  tcrossprod(get.gk.simple(X[i,], Y[i],Ti[i], 
                                            lam,be, tau1,tau0, FUNu, rho1, N, n1,...))
  }
  A<- A/N
  C<- C/N
  meat<- meat/N
  
  bread<-  matrix(0, nrow = 2*K+2, ncol = 2*K+2)
  bread[1:(2*K),1:(2*K)]<- A
  bread[2*K+1:2, ]<- cbind(C,diag(c(-1,-1)))
  bread<- solve(bread)
  #   A.inv<- solve(A)
  #   bread[-(2*K+1),-(2*K+1)]<- A.inv
  #   bread[2*K+1,]<- cbind(C%*%A.inv,-1)
  (bread%*%meat%*%t(bread))/N
}

###############################################################################
#This function obtains the big covariance matrix for all coefficients for ATT
#THe bottom right element of this will be the variance estimate for tau-ATT
get.cov.ATT<- function(X, Y, Ti, FUNu, rho, rho1, rho2, obj,...){
  N<- length(Y)
  n1<- sum(Ti)
  frac<-n1/N
  
  be<- as.vector(obj$lam.q)
  tau1<- obj$Y1
  tau0<- obj$Y0
  umat <- t(apply(X, 1, FUNu))
  K<- length(be)
  
  A<- matrix(0,ncol = K, nrow = K)
  C<- matrix(0, ncol = K, nrow = 3)
  q<- numeric(K)
  
  meat<- matrix(0, ncol = K+3, nrow = K+3)
  for(i in 1:N){
    A<- A + (1-Ti[i])*rho2(crossprod(be, umat[i,]), ...)*tcrossprod(umat[i,])
    
    C[2,(1:K)]<- C[2,(1:K)] + (1-Ti[i])*rho2(crossprod(be, umat[i,]), ...)*Y[i]*umat[i,]
    
    q <- q+1/frac^2*Ti[i]*umat[i,]
    meat<- meat +  tcrossprod(get.gk.ATT(X[i,], Y[i], Ti[i], be, tau1, tau0,
                                         FUNu, rho1, N, n1, ...))
  }
  
  A<- A/N
  C<- C/N
  meat<- meat/N
  
  bread<-  matrix(0, nrow = K+3, ncol = K+3)
  bread[1:K,1:K]<- A
  bread[K+1:3, ]<- cbind(C,diag(c(-1,-1,-1)))
  bread[1:K,K+3]<- q/N
  bread<- solve(bread)
  
  (bread%*%meat%*%t(bread))/N
}

###############################################################################
#This function obtains the BIG covariance matrix for all coefficients for MT
#This time there is no bottom right element.
get.cov.MT<- function(X, Y, Ti, FUNu, rho, rho1, rho2, obj,...){
  N<- length(Y)
  lam.mat<- obj$lam.mat
  umat <- t(apply(X, 1, FUNu))
  K<- ncol(lam.mat)
  J<- length(unique(Ti))
  taus<- obj$Yj.hat
  
  A<- matrix(0,ncol = J*K, nrow = J*K)
  C<- matrix(0, ncol = J*K, nrow = J)
  meat<- matrix(0, ncol = J*K+J, nrow = J*K+J)
  for(i in 1:N){
    for(j in 0:(J-1) ){
      temp.Ti<- 1*(Ti[i]==j)
      A[  (j*K + 1):((j + 1)*K) , (j*K + 1):((j + 1)*K) ]<- A[ (j*K + 1):((j + 1)*K) , 
                                                               (j*K + 1):((j + 1)*K) ] + 
        temp.Ti*rho2(crossprod(lam.mat[j+1,], umat[i,]), ...)*tcrossprod(umat[i,])
      C[j+1,(j*K + 1):((j + 1)*K)]<- C[j+1,(j*K + 1):((j + 1)*K)]+ 
        temp.Ti*rho2(crossprod(lam.mat[j+1,], umat[i,]), ...)*Y[i]*umat[i,]
    }
    
    meat<- meat + tcrossprod(get.gk.MT(X[i,], Y[i], Ti[i], 
                                       lam.mat[Ti[i]+1,],FUNu, rho1, J, taus,...))
  }
  
  A<- A/N
  C<- C/N
  meat<- meat/N
  
  bread1<- cbind(A, matrix(0, ncol = J, nrow = J*K))
  bread2<- cbind(C, diag(rep(-1, J)))
  bread<- rbind(bread1,bread2)
  bread<-  solve(bread)
  
  (bread%*%meat%*%t(bread))/N
}

###############################################################################
