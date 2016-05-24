dLDA <-function(xtrain,ytrain,lambda,Vinit=NULL,eps=1e-6,maxiter=1000){ 
  if (nrow(xtrain)!=length(ytrain)){
    stop("Dimensions of xtrain and ytrain don't match!")
  } 
  
  if (any(is.na(xtrain))|any(is.na(ytrain))){
    stop("Missing values are not allowed")
  }
  
  fsd=apply(xtrain,2,sd)
  if (any(fsd<1e-13)){
      stop(paste("Some features have standard deviation less than 1e-13!",sep=""))
  }
  G=max(ytrain)
  if (!is.null(Vinit)) {
      if ((nrow(Vinit)!=ncol(xtrain))|(ncol(Vinit)!=G-1)){
          stop("Supplied initial value for Vinit has wrong dimensions!")
    }
  }
  
  #center and scale X
  Xadj=scale(xtrain)
  coef=attr(Xadj,which="scaled:scale")
   
  if (G==2){
    n=length(ytrain)
    Z=matrix(0,n,2)
    for (g in 1:2){
        Z[ytrain==g,g]=1
    }
    n1=sum(ytrain==1)
    n2=n-n1
    Ytilde=sqrt(n1*n2)*Z%*%c(1/n1,-1/n2)
    V=.solveMyLasso_c(Xadj,Ytilde,lambda=lambda,eps=eps,maxiter=maxiter,binit=Vinit)
    
  } else {   
    Ytilde=.createY(ytrain)  
    V=.solveMyLassoF_c(Xadj,Ytilde,lambda=lambda,eps=eps,maxiter=maxiter,binit=Vinit) #works
  }
  
  diag(1/coef)%*%V
}
