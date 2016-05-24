KSS.CV <- function(kappa.interv, Y, X, N, T, P, spar.dim.fit, tol=tol){
  ## Y (TN x 1)
  ## X (TN x P)

  ## Estimation of Dimension d:-------------------------------------------------------------------------
  Y.mat <- matrix(Y, T,     N)					    # (T x N)
  X.mat <- matrix(X, T, (N*P))                                      # (T x N*P)
      
  Y.mat.smth  <- smooth.Pspline(x = seq(0,1,len=T), y = Y.mat, spar = spar.dim.fit)$ysmth  #(T x (N))    
  X.mat.smth  <- smooth.Pspline(x = seq(0,1,len=T), y = X.mat, spar = spar.dim.fit)$ysmth  #(T x (N)P)
           
  ## calculate beta coefficents 
  Y.smth        <- matrix(Y.mat.smth,  nrow= ((N)*T), ncol = 1)      # (T(N) x 1)
  X.smth        <- matrix(X.mat.smth,  nrow= ((N)*T), ncol = P)      # (T(N) x P)
           
  t.X.X         <- crossprod(X)       			               # (PxP)
  t.X.X.smth    <- crossprod(X, X.smth)		               # (PxP)
           
  t.X.Y         <- crossprod(X, Y)     		               # (Px1)
  t.X.Y.smth    <- crossprod(X, Y.smth)   		       # (Px1)
      
  bloc1         <- t.X.X - t.X.X.smth     		       # (PxP)
  bloc2         <- t.X.Y - t.X.Y.smth     	               # (Px1)

  ## common-Slope.Coefficients:
  com.slops.0   <- solve(bloc1)%*%bloc2				           # (Px1)
  ## calculate first step residuals and estimate dimension of factor-structure
  Residu.mat          <- matrix((Y - (X %*% com.slops.0)), T, (N)) # (Tx(N-1))        
  ## functional pca
  fpca.fit.obj     <- fpca.fit(Residu.mat, spar=spar.dim.fit)
  d.hat            <- c(OptDim(Obj=Residu.mat, criteria="KSS.C", spar=spar.dim.fit)$summary)
  ## -------------------------------------------------------------------------------------------------
  
  Outer.CV <- function(kappa){
    Inner.CV <- function(i, kappa){
      Y.mat <- matrix(Y, T,     N)					    # (T x N)
      X.mat <- matrix(X, T, (N*P))                                          # (T x N*P)
      
      Y.mat.min_i <- Y.mat[,-i]                                             # (T x N-1)
      X.mat.min_i <- X.mat[,-seq(from=c(i),to=c(N*P),by=N)]                 # (T x (N-1)*P)
      
      Y.min_i <- matrix(Y.mat.min_i, ncol=1) # (T(N-1) x 1)
      X.min_i <- matrix(X.mat.min_i, ncol=P) # (T(N-1) x P)
      
      Y.mat.min_i.smth  <- smooth.Pspline(x = seq(0,1,len=T), y = Y.mat.min_i, spar = kappa)$ysmth  #(T x (N-1))    
      X.mat.min_i.smth  <- smooth.Pspline(x = seq(0,1,len=T), y = X.mat.min_i, spar = kappa)$ysmth  #(T x (N-1)P)
           
      ## calculate beta coefficents 
      Y.min_i.smth        <- matrix(Y.mat.min_i.smth,  nrow= ((N-1)*T), ncol = 1)      # (T(N-1) x 1)
      X.min_i.smth        <- matrix(X.mat.min_i.smth,  nrow= ((N-1)*T), ncol = P)      # (T(N-1) x P)
           
      t.X.X.min_i         <- crossprod(X.min_i)       			               # (PxP)
      t.X.X.min_i.smth    <- crossprod(X.min_i, X.min_i.smth)		               # (PxP)
           
      t.X.Y.min_i         <- crossprod(X.min_i, Y.min_i)     		               # (Px1)
      t.X.Y.min_i.smth    <- crossprod(X.min_i, Y.min_i.smth)   		       # (Px1)
      
      bloc1               <- t.X.X.min_i - t.X.X.min_i.smth     		       # (PxP)
      bloc2               <- t.X.Y.min_i - t.X.Y.min_i.smth     	               # (Px1)

      ## common-Slope.Coefficients:
      com.slops.0.min_i   <- solve(bloc1)%*%bloc2				           # (Px1)
      ## calculate first step residuals and estimate dimension of factor-structure
      Residu.mat          <- matrix((Y.min_i - (X.min_i %*% com.slops.0.min_i)), T, (N-1)) # (Tx(N-1))        
      ## functional pca
      fpca.fit.obj     <- fpca.fit(Residu.mat, spar=kappa)
      Reminder_i       <- Y.mat[,i] - X.mat[,seq(from=c(i),to=c(N*P),by=N)] %*% com.slops.0.min_i      
      if(d.hat == 0){
        Sum.Resid_i    <- sum(Reminder_i^2)
      }else{
        factors        <- fpca.fit.obj$factors[,0:min(d.hat,round(sqrt(min(N, T)))), drop= FALSE]
        Sum.Resid_i    <- sum(residuals(lm(Reminder_i~factors))^2)
      }
    }
    result <-  sum(apply(matrix(1:N, N, 1),1, Inner.CV, kappa=kappa))
    cat(".")
    result
  }
  cat("Progress: CV-Optimization is running.\n")
  return.obj <- optimize(f=Outer.CV, interval=kappa.interv, tol=tol)
  return(return.obj) 
}



