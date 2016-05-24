IndivdRUMIndMH <-
function(y,X,sim=12000,burn=2000,acc=0,b0,B0,start,verbose=500){

      # start timer
      starttime=proc.time()[3]
    
      # check input    
      if(!is.vector(y)) stop("observed data y must be a vector")
      if(!is.matrix(X)) stop("the design X must be a matrix")
      if(any(is.na(y))) stop("can't handle NAs in y")
      if(any(is.na(X))) stop("can't handle NAs in X")
      if(verbose<0) stop("verbose must be a non-negative integer")
      y=as.integer(y)
      verbose=as.integer(verbose)
      
      # defining dimensions
      dims=ncol(X)                     
      N=length(y)
      if(N!=nrow(X)) stop("wrong dimensions")
      
      # parameters for the prior distribution of beta
      if(missing(b0)) b0=rep(0,dims)
      if(missing(B0)) B0=diag(10,nrow=dims,ncol=dims)
      B0inv=solve(B0)
      B0invb0=B0inv%*%b0
      
      # starting values for lambda (=exp(X*beta))
      pidach=rep(pmin(pmax(sum(y)/N,0.05),0.95),N)
      lam0=pidach/(1-pidach)
      
      # starting values for z
      U=runif(N,0,1)
      z0=log(lam0*U+y)-log(1-U+lam0*(1-y))
      
      # starting values for regression parameters beta
      if(missing(start)) start=rep(0,dims)
      
      # passing on initial values
      beta_old=beta_new=start
      lam_new=lam0                       
      z_new=z0
      
      # matrices to save the draws for beta and yi*
      beta=matrix(NA,dims,sim)
      
      # counter for acceptance rate
      count=numeric(sim)
      
      # calculation of BN for posterior distribution
      BN=solve(B0inv+3/(pi^2)*t(X)%*%X)                             
      BNinv=solve(BN)
      C=chol(BN)
      tC=t(C)
      
      # pre-calculations for bN for posterior distribution
      KN=BN%*%(3/(pi^2)*t(X))
      mN=BN%*%B0invb0
      
      # internal function: calculation of z (utilities)
      zi=function(lam_new){
             U=runif(N,0,1)
             log(lam_new*U+y)-log(1-U+lam_new*(1-y))
         }
      
      # internal function: calculation of log acceptance rate
      logprob=function(beta_new,beta_old){
                  (sum((X%*%(beta_new-beta_old)))
                   +sum(-2*log(1+exp(-z_new+X%*%beta_new)))
                   -sum(-2*log(1+exp(-z_new+X%*%beta_old)))
                   -1/2*t(beta_new-beta_old)%*%B0inv%*%(beta_new+beta_old-2*b0)
                   +1/2*t(beta_new-beta_old)%*%BNinv%*%(beta_new+beta_old-2*bN))            
              }
      
      # ----------------------- start independence MH-sampler -----------------------
      for(s in seq_len(sim)){
        
        if(s==(burn+1)){stop1=proc.time()[3]}
        
        # Step 1 - Betas
        
              bN=mN+KN%*%z_new
              beta_new=bN+tC%*%rnorm(dims,0,1)
              if(s>acc){
              u=runif(1,0,1)
              prob=exp(logprob(beta_new,beta_old))
              if(u<prob){
              count[s]=1
              } else {
              count[s]=0
              beta_new=beta_old
              }}
        
        # Step 2 - Utilities
        
              lam_new=exp(X%*%beta_new)
              z_new=zi(lam_new)
        
        # save new draws
        beta[,s]=beta_new
        beta_old=beta_new
        
        # reports
        if(verbose>0){
          if(is.element(s,c(1:5,10,20,50,100,200,500))){
            cat("sim =", s,"/ duration of iter proc so far:", 
                round(diff<-proc.time()[3]-starttime,2), "sec.,  expected time to end:", round((diff/(s-1)*sim-diff)/60,2), " min. \n")
            flush.console()
          }
          else if(s%%verbose==0){
            cat("sim =", s,"/ duration of iter proc so far:", 
                round(diff<-proc.time()[3]-starttime,2), "sec.,  expected time to end:", round((diff/(s-1)*sim-diff)/60,2), " min. \n")
            flush.console()
          }
        }
      }
      # ----------------------- end independence MH-sampler -----------------------
      
      # stop timer
      finish=proc.time()[3]
      duration=finish-starttime
      duration_wBI=finish-stop1
      
      # acceptance rate
      rate=round(sum(count[(burn+1):sim])/(sim-burn)*100,2)
      
      # output
      out <- list(beta=beta,sim=sim,burn=burn,acc=acc,dims=dims,N=N,b0=b0,B0=B0,duration=duration,
                  duration_wBI=duration_wBI,rate=rate)
      class(out) <- c("binomlogit","binomlogitMH","binomlogitIndiv")
      return(out)
}
