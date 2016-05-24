# @title Iterated Normed Gradient
# @author F. mortier, C. Trottier, G. Cornu and X. Bry
# @param Z matrix of the working variables
# @param X matrix of the normalized covariates
# @param AX matrix of the suplementary covariates
# @param W matrix of weights
# @param u vector of loadings
# @param iter_ing max number of iterations 
# @return gamma the updated loading vectors
ing <- function(Z,X,AX,W,u,method) {
  
  unew <- updateGamma(Z=Z,X=X,AX=AX,W=W,gamma=u,method=method) 
  vunew <- hFunct(z=Z,X=X,AX=AX,W=W,u=unew,method=method)$h
  crit_iterative <- 1e10
  ing_iter <- 0
  repeat{
    ustar <- updateGamma(Z=Z,X=X,AX=AX,W=W,gamma=unew,method=method)
    crit_iterative <- sum((ustar -unew)^2)
    unew <- ustar
    ing_iter <- ing_iter+1
    if(crit_iterative<method$epsilon) break
    if(ing_iter>method$maxiter) break
  }
  return(unew)
}

updateGamma <- function(Z,X,AX,W,gamma,method){
  g <- hFunct(z=Z,X=X,AX=AX,W=W,u=gamma,method=method)
  newgamma <- g$gradh/sqrt(sum(g$gradh^2))
  tmp <- hFunct(z=Z,X=X,AX=AX,W=W,u=newgamma,method=method)
  k <- 1
  while((tmp$h<g$h)&(k<method$bailout)){
    newgamma <- gamma+newgamma
    newgamma <- newgamma/sqrt(sum(newgamma^2))
    tmp <- hFunct(z=Z,X=X,AX=AX,W=W,u=newgamma,method=method)
    
    k <- k+1
  }
  return(newgamma) 
}


hFunct<- function(z,X,AX,W,u,method)
{
  #dim(W)=n*q
  #calcul des projections
  f <- c(X%*%u)
  psi <- 0
  gradpsi <- rep(0,length(u)) 
  if(!is.null(AX)){
    
    #scalsqzz <- colSums(z^2*W)
    for(k in 1:ncol(z)){
      #t(AX)%*%diag(W[,k])%*%AX
      AXtWkAX <- crossprod(AX,W[,k]*AX)  
      #Pi_T Xu = Pi_T f verif AX%*%solve(t(AX)%*%diag(W[,k])%*%AX,t(AX)%*%diag(W[,k])%*%f)
      projWkfAX <- c(AX%*%solve(AXtWkAX,crossprod(AX,W[,k]*f)))
      # Pi_{T^ortho} Xu = Pi_{T^ortho} f 
      projWkforthoAX <- f - projWkfAX  
      #Pi_T z_k verif AX%*%solve(t(AX)%*%diag(W[,k])%*%AX,t(AX)%*%diag(W[,k])%*%z[,k])
      #z W_k standardized
      zk <- wtScale(z[,k],W[,k])#z[,k] - sum(W[,k]*z[,k])
      Wzk <- W[,k]*zk
      projWkzAX <- AX%*%solve(AXtWkAX,crossprod(AX,Wzk))      
      #projWkzorthoAX <- z[,k] - projWkzAX #Pi_{T^{ortho}} z_k      
      #projWkzorthoAX <- Wz[,k] - projWkzAX  
      #calcul de psi
      scalsqpfz <- sum(c(projWkforthoAX)*Wzk)^2#<Pi_{T^{ortho}}Xu|z_k>^2_{W_k} verif  (t(projWkforthoAX)%*%diag(W[,k])%*%z[,k])^2
      scalsqpfpf <- sum(c(projWkforthoAX)^2*W[,k])#||Pi_{T^{ortho}}Xu||_{W_k}^2    
      #term1psi <- sum(scalsqpfz/(scalsqpfpf*scalsqzz[k]))
      ##comme z[,k] est W[,k] standardized, le terme scalsqzz[k] disparait
      term1psi <- sum(scalsqpfz/(scalsqpfpf))
      term2psi <- sum(Wzk*projWkzAX)
      psi <- psi+term1psi+term2psi     
      #calcul de grad de psi
      PiorthoPrimeWkz <- Wzk -  W[,k]*AX%*%solve(AXtWkAX,crossprod(AX,Wzk))
      ##verif Wz[,k] - diag(W[,k])%*%AX%*%solve(t(AX)%*%diag(W[,k])%*%AX)%*%t(AX)%*%diag(W[,k])%*%z[,k]
      XprimeprojorthoWz <- crossprod(X,PiorthoPrimeWkz) #X' Pi_{T^{ortho}}'W_k z_k
      term1 <- c(XprimeprojorthoWz%*%crossprod(XprimeprojorthoWz,u))/(scalsqpfpf)
      
      WprojWkOrthof <- W[,k]*projWkforthoAX#W_k Pi_{T^{ortho}}Xu
      ##verif diag(W[,k])%*%projWkforthoAX
      # PiorthoPrimeWkf <- WprojWkOrthof-W[,k]*AX%*%solve(AXtWkAX,crossprod(AX,WprojWkOrthof))#Pi_{T^{ortho}}^primeW_k\pi_{T^{ortho}}Xu
      
      ##verif cf  PiorthoPrimeWkz
      #term2 <-  scalsqpfz*c(crossprod(X,PiorthoPrimeWkf))/(scalsqpfpf^2*scalsqzz[k])  
      
      term2 <-  scalsqpfz*c(crossprod(X,WprojWkOrthof))/(scalsqpfpf^2)  
      gradpsi <- gradpsi +(term1-term2)
    }
    gradpsi <- 2*gradpsi
  }else{
    for(k in 1:ncol(z)){
      zk <- wtScale(z[,k],W[,k])#z[,k] - sum(W[,k]*z[,k])
      Wzk <- W[,k]*zk
      scalsqpfz <- sum(c(f)*Wzk)^2
      scalsqpfpf <- sum(c(f)^2*W[,k])
      #calcul de psi
      psi <- psi+sum(scalsqpfz/(scalsqpfpf))      
      #calcul de grad de psi     
      XprimeWz <- crossprod(X,Wzk) #X'W_k z_k
      term1 <- c(XprimeWz%*%crossprod(XprimeWz,u))/(scalsqpfpf)
      term2 <-  scalsqpfz*c(crossprod(X,W[,k]*f))/(scalsqpfpf^2)  
      gradpsi <- gradpsi +(term1-term2)
    }
    gradpsi <- 2*gradpsi
  }
  n <- nrow(X)    
  # calcul phi Component Variance: cv
  if(method$phi=="cv") {
    phi <- c(crossprod(f))/n    
    # calcul grad phi
    gradphi <- c(2*crossprod(X,f/n))
  } else { 
    ### autre calcul de phi avec l>=1 : vpi: Variable Powered Inertia
    scalsqfX <- colSums(f*X/n)
    XtWX <- crossprod(X)/n
    phi <- (sum((scalsqfX^2)^method$l))^(1/method$l)    
    # calcul de grad phi
    gradphi <- 2*phi^(1-method$l)*rowSums(XtWX%*%diag(scalsqfX)^(2*method$l-1))
  }
  # calcul de h (s in R+)
  #h = log(psi)+method$s*log(phi)
  #gradh=gradpsi/psi+method$s*gradphi/phi
  # calcul de h (s in [0..1])
  h = (1-method$s)*log(psi)+method$s*log(phi)
  gradh=(1-method$s)*gradpsi/psi+method$s*gradphi/phi    
  return(list(h=h, gradh=gradh,psi=psi,gradpsi=gradpsi,phi=phi,gradphi=gradphi))
}


