krls <-
function(     X=NULL,
              y=NULL,
              whichkernel="gaussian",
              lambda=NULL,
              sigma=NULL,
              derivative=TRUE,
              binary=TRUE,
              vcov=TRUE,
              print.level=1,
              L=NULL,
              U=NULL,
              tol=NULL,
							eigtrunc=NULL){
              
      # checks
      y <- as.matrix(y)
      X <- as.matrix(X)
      
      if (is.numeric(X)==FALSE){
       stop("X must be numeric")
      }      
      if (is.numeric(y)==FALSE){
       stop("y must be numeric")
      }
      if (sd(y)==0) {
        stop("y is a constant")        
      }  
       if (sum(is.na(X))>0){
       stop("X contains missing data")
      }
      if (sum(is.na(y))>0){
       stop("y contains missing data")
      }
      if (var(y)==0){
       stop("y does not vary")
      }
      
      if (!is.null(eigtrunc)){
          if (!is.numeric(eigtrunc)) stop("eigtrunc, if used, must be numeric")    
          if (eigtrunc>1 | eigtrunc<0) stop("eigtrunc must be between 0 and 1")
          if (eigtrunc==0) {
            eigtrunc=NULL
            warning("eigtrunc of 0 equivalent to no eigen truncation")}
      }
      
      n <- nrow(X)
      d <- ncol(X)
      if (n!=nrow(y)){
       stop("nrow(X) not equal to number of elements in y")
      }
      
      stopifnot(
                is.logical(derivative),
                is.logical(vcov),
                is.logical(binary)
                )
      
      if(derivative==TRUE){
        if(vcov==FALSE){
        stop("derivative==TRUE requires vcov=TRUE")
        }
      }
      
      # default sigma to dim of X 
      if(is.null(sigma)) { sigma <- d
      } else {
        stopifnot(is.vector(sigma),
                  length(sigma)==1,
                  is.numeric(sigma),
                  sigma>0)        
      }
      
      # column names
      if (is.null(colnames(X))) {
      colnames(X) <- paste("x",1:d,sep="")
      }
         
      # scale
      X.init <- X
      X.init.sd <- apply(X.init,2,sd)
      if (sum(X.init.sd==0)) {
        stop("at least one column in X is a constant, please remove the constant(s)")
      }      
      y.init <- y
      y.init.sd <- apply(y.init,2,sd)
      y.init.mean <- mean(y.init)
      X <- scale(X,center=TRUE,scale=X.init.sd)    
      y <- scale(y,center=y.init.mean,scale=y.init.sd)
     
      # kernel matrix
      K <- NULL
      if(whichkernel=="gaussian"){ K <- gausskernel(X,sigma=sigma)}
      if(whichkernel=="linear"){K <- tcrossprod(X)}
      if(whichkernel=="poly2"){K <- (tcrossprod(X)+1)^2}
      if(whichkernel=="poly3"){K <- (tcrossprod(X)+1)^3}
      if(whichkernel=="poly4"){K <- (tcrossprod(X)+1)^4}
      if(is.null(K)){stop("No valid Kernel specified")}
      
      # eigenvalue decomposition
      Eigenobject=list()
      Eigenobject$values=numeric(length=n)
      Eigenobject$vectors=matrix(NA,n,n)
      Eigenobject <- eigen(K,symmetric=TRUE)
      
      # default lamda is chosen by leave one out optimization 
       if(is.null(lambda)) {
 
 #  Old Way to get lambda       
 # first try with max eigenvalue (increase interval in case of corner solution at upper bound)
 #   lowerb  <-  .Machine$double.eps
 #   upperb <- max(Eigenobject$values)
 #     if(print.level>1) { cat("Using Leave one out validation to determine lamnda. Search Interval: 0 to",round(upperb,3), "\n")}  
 #     lambda <- optimize(looloss,interval=c(lowerb,upperb),y=y,Eigenobject=Eigenobject)$minimum
 #     if(lambda >= (upperb - .5)){
 #       if(print.level>1) { cat("Increasing search window for Lambda that minimizes Loo-Loss \n")}  
 #      lambda <- optimize(looloss,interval=c(upperb,2*upperb),y=y,Eigenobject=Eigenobject)$minimum
 #      }
      # run golden section search
       if(print.level>2){noisy <- TRUE
       } else {
         noisy <- FALSE
       }   
       lambda<- lambdasearch(L=L,U=U,y=y,Eigenobject=Eigenobject,eigtrunc=eigtrunc,noisy=noisy)
         
       if(print.level>1) { cat("Lambda that minimizes Loo-Loss is:",round(lambda,5),"\n")}    
       
       } else {  # check user specified lambda
         stopifnot(is.vector(lambda),
                   length(lambda)==1,
                   is.numeric(lambda),
                   lambda>0)  
      }
      # solve given LOO optimal or user specified lambda
      out <-  solveforc(y=y,Eigenobject=Eigenobject,lambda=lambda,eigtrunc=eigtrunc)
      
 #Now, if eigtrunc was used, maybe we should also reconstruct K using those eigenvalues?
 #if (!is.null(eigtrunc)){	
 #lastkeeper=max(which(Eigenobject$values>=eigtrunc*Eigenobject$values[1]))	
 #K <- tcrossprod(multdiag(X=Eigenobject$vectors[,1:lastkeeper],d=Eigenobject$values[1:lastkeeper]),Eigenobject$vectors[,1:lastkeeper])        
 #}

 # fitted values
      yfitted <- K%*%out$coeffs

      ## var-covar matrix for c
      if(vcov==TRUE){      
        # sigma squared
        sigmasq <- as.vector((1/n) * crossprod(y-yfitted))
        #Gdiag   <-  Eigenobject$values+lambda
        #Ginvsq  <-  tcrossprod(Eigenobject$vectors %*% diag(Gdiag^-2,n,n),Eigenobject$vectors)  
        #Faster:
        #Ginvsq <- tcrossprod(multdiag(X=Eigenobject$vectors,d=Gdiag^-2),Eigenobject$vectors)

        # Faster still: skip creating Ginvsq and multiply by sigmasq at same time
        vcovmatc <- matrix(NA,n,n)
        
        if (is.null(eigtrunc)){	
        	vcovmatc <- tcrossprod(multdiag(X=Eigenobject$vectors,d=sigmasq*(Eigenobject$values+lambda)^-2),Eigenobject$vectors)        
        } else{
        	
        	#eigentruncation: keep only eigenvectors at least 'eigtrunc' times as large as the largest
        	lastkeeper=max(which(Eigenobject$values>=eigtrunc*Eigenobject$values[1]))	
        	vcovmatc <- tcrossprod(multdiag(X=Eigenobject$vectors[,1:lastkeeper],d=sigmasq*(Eigenobject$values[1:lastkeeper]+lambda)^-2),Eigenobject$vectors[,1:lastkeeper])        
        }
        
        #vcovmatc <- sigmasq*Ginvsq
        
				# var-covar for y hats
				vcovmatyhat <- crossprod(K,vcovmatc%*%K) 

      } else { 
        vcov.c      <- NULL
        vcov.fitted <- NULL
      }
                
      # compute derivatives
      avgderiv <- varavgderivmat <- derivmat <- NULL
   
 			if(derivative==TRUE){
        if(whichkernel!="gaussian"){
          stop("derivatives are only available when whichkernel='gaussian' is specified")
        }
        
      derivmat<-matrix(NA,n,d)
      varavgderivmat<- avgderivmat <- matrix(NA,1,d)
        
      rows <- cbind(rep(1:nrow(X), each = nrow(X)), 1:nrow(X))
      distances <- X[rows[,1],] - X[ rows[,2],]    # d by n*n matrix of pairwise distances  
      colnames(derivmat)       <- colnames(X)
      colnames(varavgderivmat) <- colnames(X)
      
      for(k in 1:d){       
       if(d==1){
             distk <-  matrix(distances,n,n,byrow=TRUE)
         } else {
             distk <-  matrix(distances[,k],n,n,byrow=TRUE) 
        }
         L <-  distk*K
         # pointwise derivatives
         derivmat[,k] <- (-2/sigma)*L%*%out$coeff
         # variance for average derivative
         varavgderivmat[1,k] <- (1/n^2)*sum((-2/sigma)^2 * crossprod(L,vcovmatc%*%L))     
      }
       # avg derivatives
       avgderiv <- matrix(colMeans(derivmat),nrow=1)
       colnames(avgderiv) <- colnames(X)
      # get back to scale
      derivmat <- scale(y.init.sd*derivmat,center=FALSE,scale=X.init.sd)
      attr(derivmat,"scaled:scale")<- NULL
      avgderiv <- scale(as.matrix(y.init.sd*avgderiv),center=FALSE,scale=X.init.sd)
      attr(avgderiv,"scaled:scale")<- NULL
      varavgderivmat <- (y.init.sd/X.init.sd)^2*varavgderivmat
      attr(varavgderivmat,"scaled:scale")<- NULL
      }
         
      # get back to scale
      yfitted     <- yfitted*y.init.sd+y.init.mean
       if(vcov==TRUE){  
           vcov.c      <- (y.init.sd^2)*vcovmatc
           vcov.fitted <- (y.init.sd^2)*vcovmatyhat
       } else {
         vcov.c      <- NULL
         vcov.fitted <- NULL
      }
      Looe        <- out$Le*y.init.sd
      # R square
      R2 <- 1-(var(y.init-yfitted)/(y.init.sd^2))
            
      # indicator for binary predictors
      binaryindicator=matrix(FALSE,1,d)
      colnames(binaryindicator) <- colnames(X)
      
      # return
   z <- list(K=K,
             coeffs=out$coeffs,
             Looe=Looe,
             fitted=yfitted,
             X=X.init,
             y=y.init,
             sigma=sigma,
             lambda=lambda,
             R2 = R2,
             derivatives=derivmat,
             avgderivatives=avgderiv, 
             var.avgderivatives=varavgderivmat,
             vcov.c=vcov.c,
             vcov.fitted=vcov.fitted,
             binaryindicator=binaryindicator
            )          
  class(z) <- "krls"
    
  # add first differences if requested
    if(derivative==TRUE && binary==TRUE){
      z <- fdskrls(z)
    }
  
  # printing
      if(print.level>0 && derivative==TRUE){
      output <- setNames(as.vector(z$avgderivatives), colnames(z$avgderivatives))
      cat("\n Average Marginal Effects:\n \n")
      print(output)
      
      cat("\n Quartiles of Marginal Effects:\n \n")
      print(apply(z$derivatives,2,quantile,probs=c(.25,.5,.75)))
      }
      
  return(z)
            
}

