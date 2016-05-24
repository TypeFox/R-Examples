# Y: A metabolomics data matrix with samples in rows and metabolites in columns
# ctl: A logical vector indicating the controls
# lambda: regularization parameter for the unwanted variation.
# k: rank of the unwanted variation term.

RUVRand <- function(Y, ctl,lambda=NULL, k=NULL,plotk=TRUE,...){
  
  Yc<-Y[, ctl]
  svdYc <- svd(Yc)
  fullW <- svdYc$u %*% diag(svdYc$d)  
  if(plotk){
    barplot(prcomp(Yc,scale. =T)$sdev^2/
              sum(prcomp(Yc,scale. =T)$sdev^2),
#            xlim=c(0,10),
#            names.arg =c(1:9),
            ylim=c(0,1),
            xlab="k",
            ylab="proportion of the variance",
            cex.lab=1.2,cex.axis=1.2,las=2)
    
  }
  if(is.null(k))
    stop('k must be entered')    
  if (!is.null(k) & is.null(lambda)){ 
    optklambda<-opt(ktry=k, W=fullW,Yc=Yc)
    lambda<-optklambda$optmat[,3]     
  } else optklambda<-NULL
  
  W<-fullW[,1:k,drop=FALSE]
  alpha<-solve(t(W)%*%W + lambda*diag(k), t(W) %*% Y)
  uvcomp<-W %*% alpha
  newY <- Y - uvcomp 
return(list(unadjY=Y,newY=newY,UVcomp=uvcomp,W=W,alpha= alpha,opt=optklambda,
            k=k,lambda=lambda,ctl=ctl))  
}


# Y: A metabolomics data matrix with samples in rows and metabolites in columns
# W: A data matrix of unwanted variation with samples in rows and factors in columns
loglik<- function (par, Y,W) 
{
  m <- ncol(Y)
  n<-nrow(Y)
  if (!is.null(W)){
    sigma2.a<-par[1]
    sigma2.e<-par[2]  
    if ((sigma2.a<0)|(sigma2.e<0))
      return(1e6)
    Sigma<-sigma2.a*(W%*%t(W))+sigma2.e*diag(m)
  } else{
    Sigma <- diag(m) 
    Sigma[upper.tri(Sigma, diag=TRUE)] <- par 
    Sigma <- Sigma + t(Sigma) - diag(diag(Sigma)) 
  }
  ed = eigen(Sigma, symmetric = TRUE)
  ev = ed$values
  if (!all(ev >= -1e-06 * abs(ev[1])))
    return(1e6)
  mu<-rep(0,m)
  centeredx<-sweep(Y,2,mu,"-")
  ssnew<-  t(centeredx)%*%(centeredx)
  if (is.null(tryCatch(solve(Sigma), error=function(e) NULL)))
    return(1e6)
  else
    inv.Sigma<-solve(Sigma) 
  Sigmainvss<-inv.Sigma%*%ssnew
  return(n*determinant(Sigma,logarithm=T)$mod+sum(diag(Sigmainvss)))
  
}

#kvec : A numerical vector with values of k for which lambda needs to be estimated
# W: A data matrix of unwanted variation with samples in rows and factors in columns
# Yc: A data matrix with samples in rows and quality control metabolites in columns
opt<-function(ktry,W,Yc){
  opt<-list()
  optmat<-matrix(NA,nrow=1,ncol=8)
  colnames(optmat)<-c("sigma2.a","sigma2.e","nu",
                      "lower_sigma2.a","upper_sigma2.a",
                      "lower_sigma2.e","upper_sigma2.e",
                      "convergence")
  
    opt<-optim(c(0.1,0.1),
                    loglik,
                    Y=t(Yc),
                    W=W[,1:ktry,drop=FALSE],
                    hessian=T)    
    fisher_info<-solve(opt$hessian/2)
    se<-sqrt(diag(fisher_info))
    upper_par1<-opt$par[1]+1.96*se[1]
    lower_par1<-opt$par[1]-1.96*se[1]
    upper_par2<-opt$par[2]+1.96*se[2]
    lower_par2<-opt$par[2]-1.96*se[2]
    
    optmat[1,]<-c(opt$par[1],
                  opt$par[2],
                  opt$par[2]/opt$par[1],
                  lower_par1, upper_par1,
                  lower_par2, upper_par2,                  
                  opt$convergence)
  
  rownames(optmat)<-ktry
  return(list(optmat=optmat, opt=opt))
}

#RUVRand: Output from RUVRand
#maxIter: Maximum no of iterations
#wUpdate: Update W every wUpdate iterations
#... :Other parameters for kmeans
RuvRandIter <- function(RUVRand,
                        maxIter, 
                        wUpdate=maxIter+1, 
                        lambdaUpdate=TRUE,
                        p=p,...){
  Y<-RUVRand$unadjY
  ctl<-RUVRand$ctl
  ndim<-RUVRand$ndim
  m <- nrow(Y) 
  n <- ncol(Y) 
  
  converged <- 0
  iter <- 0
  currObj <- Inf
  cEps<-1e-6
  

  W <- RUVRand$W
  a <- RUVRand$alpha
  lambda<-RUVRand$lambda
  k<-RUVRand$k
  Wa <- W %*% a
  
  X <- matrix(0,m,p)
  b <- matrix(0,p,n)
  Xb <- X %*% b
  
  while(!converged){
    iter <- iter + 1
    
    print('Estimating the factor of interest')
    XbOld <- Xb
    
    kmres <- kmeans((Y[, -ctl] - Wa[, -ctl]),centers=p,nstart=20,...)
    idx <- kmres$cluster
    for(kk in 1:p){
      X[, kk] <- cbind(as.numeric(idx==kk))
    }
    b[, -ctl] <- kmres$centers
    Xb <- X %*% b
     
    WaOld <- Wa
    WOld <- W    
    
    if(iter / wUpdate == iter %/% wUpdate){
      
      print('Re-estimating W')
      
      svdYmXb <- svd((Y - Xb))
      fullW <- svdYmXb$u %*% diag(svdYmXb$d) 
      if (lambdaUpdate){
        print('Re-estimating k and lambda')
      barplot(prcomp((Y - Xb),scale. =T)$sdev^2/
                sum(prcomp((Y - Xb),scale. =T)$sdev^2),
              xlim=c(0,min(dim(Y))+1),
              names.arg =c(1:min(dim(Y))),
              ylim=c(0,1),
              xlab="k",
              ylab="proportion of the variance",
              cex.lab=1.2,cex.axis=1.2)
      
      k<-readk()
      W <- fullW[,c(1:k)]    
      optklambda<-opt(ktry=k,
                     W=W,
                     Yc=(Y - Xb))
      lambda<-optklambda$optmat[,3] 
      print(paste("lambda =" ,lambda))      
      }
    }
    
    print('Estimating the unwanted variation component')
    
    
    
    a <- solve(t(W)%*%W + lambda*diag(ncol(W)), t(W) %*% (Y - Xb))
    Wa <- W %*% a
    
    print('Update done')
    
  
    l2Err <- (norm((Y - Xb - Wa), 'F')^2)/(m*n)
    oldObj <- currObj
    currObj <- norm((Y - Xb - Wa), 'F')^2 +lambda*norm(a,'F')^2
    
    dXb <- norm(Xb-XbOld,'F')/norm(Xb,'F')
    dWa <- norm(Wa-WaOld,'F')/norm(Wa,'F')
    if(ncol(W) != ncol(WOld)){
      dW <- Inf}
    else{
      dW <- norm(W-WOld,'F')/norm(W,'F')      
    }  
    dObj <- (oldObj-currObj)/oldObj
    
    print(sprintf('iter %d, dXb/Xb=%g, dWa/Wa=%g, dW/W=%g,l2Err=%g, dObj/obj=%g, obj=%g',
                  iter,dXb,dWa,dW,l2Err,dObj,currObj))
    
    
    if(iter >= maxIter || (!is.nan(max(dXb,dWa)) && max(dXb,dWa) < cEps)){
      converged = 1
    }
  }
  
  cY <- Y - Wa
  
  return(list(unadjY=RUVRand$unadjY, newY=cY,  UVcomp=Wa,
              W=W,alpha=a,X=X,b=b,k=k,lambda=lambda,opt=RUVRand$opt, 
              ctl=RUVRand$ctl))
}

readk <- function()
{ 
  n <- readline(prompt="Enter k: ")
  if(!grepl("^[0-9]+$",n))
  {
    return(readk())
  }
  
  return(as.integer(n))
}
