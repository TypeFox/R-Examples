# @title Function that calculates the current scglr-component
# @export 
# @param X a matrix (n*p) containing covariates standardized
# @param Y a matrix (n*q) containing dependent variables
# @param AX a matrix containing additional covariates 
# @param family a vector of charater of length q specifying the distributions of the responses. 
# Bernoulli, binomial, Poisson and gaussian are allowed.
# @param size specifies the numbers of trials of the binomial responses.  A (n*qb) matrix is expected
#  for qb binomial variables.
# @param offset used for the poisson dependent variables.
# A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected
# @param ds an integer specifying the degree to which structural strength in X is taken into account
# when calculating the component
# @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and 
# the tolerance convergence criterion for the Fisher scoring algorithm. Default is set to 50 and 10e-6 respectively. 
# @param method optimization algorithm used for loadings and components. Object of class "method.SCGLR" 
# built by \code{\link{methodEigen}} or \code{\link{methodIng}}
# @return the unit vector of loadings associated with the current component, 
# i.e. the coefficients of the regressors in the linear combination giving each component
oneComponent <- function(X,Y,AX,family,size=NULL,offset=NULL,ds,crit,method)
{
  ##cst control and iteration
  muinf <- 1e-5
  maxit <- crit$maxit
  tol  <- crit$tol
  tol1<-1
  iter<-0
  ##dimension 
  n<-dim(X)[1]
  p<-dim(X)[2]
  q<-dim(Y)[2]
  if (is.null(q)) q<-1
  # initialisation de Z et de W, de alpha, de u et de a
  if(!is.null(offset)) loffset <- log(offset)  
  ###Initialization, Z working variables  
  mu0 <- apply(Y,2,mean)   
  mu0 <- matrix(mu0,n,q,byrow=TRUE)
  Z <- Y
  
  if("bernoulli"%in%family)
  {
    tmu0 <- mu0[,family=="bernoulli"]
    Z[,family=="bernoulli"] <-  log(tmu0/(1-tmu0))+(Y[,family=="bernoulli"]-tmu0)/(tmu0*(1-tmu0))
  }  
  if("binomial"%in%family)
  {
    tmu0 <- mu0[,family=="binomial"]
    Z[,family=="binomial"] <-   log(tmu0/(1-tmu0))+(Y[,family=="binomial"]-tmu0)/(tmu0*(1-tmu0))
  }
  if("poisson"%in%family)
  {
    tmu0 <- mu0[,family=="poisson"]
    if(is.null(offset)){
      Z[,family=="poisson"]<- log(tmu0)+(Y[,family=="poisson"]-tmu0)/tmu0
    }else{
      Z[,family=="poisson"]<- log(tmu0)-loffset+(Y[,family=="poisson"]-tmu0)/tmu0
    }
  }
  if("gaussian"%in%family){
    Z[,family=="gaussian"] <- Y[,family=="gaussian"]
  }
  
  
  Zjcr <- apply(Z,2,function(x) wtScale(x,1/n))
  if(p>n){
    SA <- tcrossprod(Zjcr)
    A <-  crossprod(X,SA)%*%X
  }else{
    SA <- crossprod(Zjcr,X)
    A <- crossprod(SA)
  }
  
  if(ds>0)
  {
    stX <- crossprod(X)/n
    puis.stX <- stX%^%ds
    A <- puis.stX%*%A
  }
  # calcul du premier u	
  u<-eigen(A,symmetric=TRUE)$vectors[,1]
  
  #if (is.logical(u)) {return(FALSE)}  
  #calcul de f, 
  f<-X%*%u 
  # initialisation des eta
  if (is.null(AX))
  {
    reg<-cbind(1,f)
    sol <-solve(crossprod(reg),crossprod(reg,Z)) 
    eta <- apply(sol,2,function(x) x[1]+X%*%(u*x[2]))				
  }
  else{
    r <- dim(AX)[2]
    reg <- cbind(1,AX,f)
    sol <-solve(crossprod(reg),crossprod(reg,Z))
    eta <- apply(sol,2,function(x) x[1]+AX%*%x[2:(r+1)]+X%*%(u*x[r+2]))
  }
  # Update initialization of Z and initialization of W 
  #Z <- eta  
  W <- matrix(1,n,q)
  
  if("bernoulli"%in%family)
  {
    etainf <- log(muinf/(1-muinf))
    indinf<-1*(eta[,family=="bernoulli"]<etainf)
    eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indinf)+etainf*indinf
    indsup<-1*(eta[,family=="bernoulli"]>-etainf)
    eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indsup)-etainf*indsup
    mu <- exp(eta[,family=="bernoulli"])/(1+exp(eta[,family=="bernoulli"]))
    Z[,family=="bernoulli"] <-  eta[,family=="bernoulli"] + (Y[,family=="bernoulli"]-mu)/(mu*(1-mu))
    W[,family=="bernoulli"] <- mu*(1-mu)
  }
  if("binomial"%in%family)
  {
    etainf <- log(muinf/(1-muinf))
    indinf<-1*(eta[,family=="binomial"]<etainf)
    eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indinf)+etainf*indinf
    indsup<-1*(eta[,family=="binomial"]>-etainf)
    eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indsup)-etainf*indsup
    mu <- exp(eta[,family=="binomial"])/(1+exp(eta[,family=="binomial"]))
    Z[,family=="binomial"] <-  eta[,family=="binomial"] + (Y[,family=="binomial"]-mu)/(mu*(1-mu))
    W[,family=="binomial"] <- mu*(1-mu)*size
  }
  if("poisson"%in%family)
  {
    etainf <- log(muinf)
    indinf<-1*(eta[,family=="poisson"]<etainf)
    eta[,family=="poisson"]<-eta[,family=="poisson"]*(1-indinf)+etainf*indinf
    if(is.null(offset))
    {
      mu <- exp(eta[,family=="poisson"])
    }else{
      mu <- exp(eta[,family=="poisson"]+loffset)   
    }
    Z[,family=="poisson"]<- eta[,family=="poisson"]+(Y[,family=="poisson"]-mu)/mu
    W[,family=="poisson"] <- mu
  } 
  if("gaussian"%in%family){
    Z[,family=="gaussian"]<-Y[,family=="gaussian"]
  }
  W <- apply(W,2,function(x) x/sum(x)) 
  
  ###Loop 
  
  if(FALSE) { ## TODO verbose option
    message(ds," ",appendLF=FALSE)
    progress <- txtProgressBar(0,maxit,style=1)
  }
  
  while((tol1>tol)&&(iter<maxit))
  {
    if(method$method=="lpls") {
      ##################################### LPLS METHOD #############
      # centrage des Zj par rapport a Wj
      # centrage de X par rapport a Wbar
      # calcul de A : matrice de la forme quadratique du critere a miximiser
      A<-matrix(0,p,p)
      # calcul de A
      if (is.null(AX)){
        for (j in 1:q){
          Xjc <- apply(X,2, wtCenter,w=W[,j]) # seulement centre ! 
          Zjcr <- wtScale(Z[,j],W[,j]) # centre-reduit !
          WZj <- W[,j]*Zjcr
          if(p>n){
            SA <- tcrossprod(WZj,WZj)
            A <- A+crossprod(Xjc,SA)%*%Xjc
          }else{
            SA <- crossprod(WZj,Xjc)
            A <- A+ crossprod(SA)
          }
        }
      }
      else{
        for (j in 1:q){
          Xjc <- apply(X,2,wtCenter,w=W[,j]) ## BOTTLENECK
          Xjc <- Xjc-AX%*%solve(crossprod(AX,W[,j]*AX),crossprod(AX,W[,j]*Xjc))
          Zjcr <- wtScale(Z[,j],W[,j])
          WZj <- W[,j]*Zjcr    
          if(p>n){
            SA <- tcrossprod(WZj,WZj)
            A <- A+crossprod(Xjc,SA)%*%Xjc
          }else{
            SA <- crossprod(WZj,Xjc)
            A <- A+ crossprod(SA)
          }        
        }
        
      }
      
      if(ds>0){
        A <- puis.stX%*%A
      }
      unew <- eigen(A,symmetric=TRUE)$vectors[,1]
      if (c(crossprod(unew,u)) < 0) {
        unew <- -unew
      } 
    } else {
      ############################ SR METHOD ########################      
      unew <- ing(Z=Z,X=X,AX=AX,W=W,u=u,method=method) 
    }
    f<-X%*%unew
    if (is.null(AX)) {
      reg <- cbind(1,f)      
      for(j in 1:q) {	
        sol <-solve(crossprod(reg,W[,j]*reg),crossprod(reg,W[,j]*Z[,j]))
        eta[,j] <- sol[1] + f%*%sol[2] 
      }		
    } else {	
      r<-dim(AX)[2]
      reg<-cbind(rep(1,n),AX,f)
      for (j in 1:q) {       
        sol <-solve(crossprod(reg,W[,j]*reg),crossprod(reg,W[,j]*Z[,j]))
        eta[,j] <- sol[1] + AX%*%sol[2:(r+1)] + f%*%sol[r+2] 
      }
    }
    
    
    # Update  of Z and W 
    #browser()
    #Z <- eta
    if("bernoulli" %in% family) {
      etainf <- log(muinf/(1-muinf))
      indinf<-1*(eta[,family=="bernoulli"]<etainf)
      eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indinf)+etainf*indinf
      indsup<-1*(eta[,family=="bernoulli"]>-etainf)
      eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indsup)-etainf*indsup
      mu <- exp(eta[,family=="bernoulli"])/(1+exp(eta[,family=="bernoulli"]))
      Z[,family=="bernoulli"] <-  eta[,family=="bernoulli"] + (Y[,family=="bernoulli"]-mu)/(mu*(1-mu))
      W[,family=="bernoulli"] <- mu*(1-mu)
    }
    if("binomial" %in% family) {
      etainf <- log(muinf/(1-muinf))
      indinf<-1*(eta[,family=="binomial"]<etainf)
      eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indinf)+etainf*indinf
      indsup<-1*(eta[,family=="binomial"]>-etainf)
      eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indsup)-etainf*indsup
      mu <- exp(eta[,family=="binomial"])/(1+exp(eta[,family=="binomial"]))
      Z[,family=="binomial"] <-  eta[,family=="binomial"] + (Y[,family=="binomial"]-mu)/(mu*(1-mu))
      W[,family=="binomial"] <- mu*(1-mu)*size
    }
    if("poisson" %in% family) {
      etainf <- log(muinf)
      indinf<-1*(eta[,family=="poisson"]<etainf)
      eta[,family=="poisson"]<-eta[,family=="poisson"]*(1-indinf)+etainf*indinf
      if(is.null(offset)) {
        mu <- exp(eta[,family=="poisson"])
      } else {
        mu <- exp(eta[,family=="poisson"]+loffset)   
      }
      Z[,family=="poisson"]<- eta[,family=="poisson"]+(Y[,family=="poisson"]-mu)/mu
      W[,family=="poisson"] <- mu
    }
    if("gaussian"%in%family) {
      Z[,family=="gaussian"]<-Y[,family=="gaussian"]
    }
    
    W <- apply(W,2,function(x) x/sum(x))   
    
    f<-c(X%*%u)
    fnew<-c(X%*%unew)
    f<-f/sqrt(c(crossprod(f)))
    fnew<-fnew/sqrt(c(crossprod(fnew)))
    tol1<-1-crossprod(f,fnew)^2 
    # TODO verbose option
    if(FALSE) setTxtProgressBar(progress, iter)
    u<-unew
    iter<-iter+1
  }
  # TODO verbose option
  if(FALSE)
    close(progress)
  
  if(iter==maxit)
  {
    # TODO verbose option
    if(FALSE)
      message("  Warning !!! max number of iterations in oneComponent tol=", tol1)
    return(FALSE)
  }
  if(tol1 <= tol) {
    # TODO verbose option
    if(FALSE) 
      message("  Convergence in ",iter," iterations with ds=", ds, " and tol=", tol)    
  }
  return(u=u)
}
