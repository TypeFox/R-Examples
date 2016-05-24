dRUMIndMH <-
function(yi,Ni,X,sim=12000,burn=2000,acc=0,b0,B0,start,verbose=500){
  
      # start timer
      starttime=proc.time()[3]
      
      # check input
      if(!is.vector(yi)) stop("the counts yi must be a vector")
      if(!is.vector(Ni)) stop("the numbers of trials Ni must be a vector")
      if(!is.matrix(X)) stop("the design X must be a matrix")
      if(length(yi)!=length(Ni)) stop("yi and Ni must have same length")
      if(any(is.na(yi))) stop("can't handle NAs in yi")
      if(any(is.na(Ni))) stop("can't handle NAs in Ni")
      if(any(is.na(X))) stop("can't handle NAs in X")
      if(verbose<0) stop("verbose must be a non-negative integer")
      yi=as.integer(yi)
      Ni=as.integer(Ni)
      verbose=as.integer(verbose)
      
      # delete lines, when covariable pattern was not observed (i.e. Ni=0)
      if(any(Ni==0)){
        X=X[(Ni>0),,drop=FALSE]
        yi=yi[(Ni>0)]
        Ni=Ni[(Ni>0)]
      }
      
      # defining dimensions
      dims=ncol(X)                     
      t=length(yi)
      if(t!=nrow(X)) stop("wrong dimensions") 

      # parameters for the prior distribution of beta
      if(missing(b0)) b0=rep(0,dims)
      if(missing(B0)) B0=diag(10,nrow=dims,ncol=dims)
      B0inv=solve(B0)
      B0invb0=B0inv%*%b0                                  
      
      # indicators used for calculation of yi*
      indi1=as.numeric(yi>0)
      indi2=as.numeric(yi<Ni)

      # starting values for lambda (=exp(X*beta))
      pidach=pmin(pmax(yi/Ni,0.05),0.95)
      lam0=pidach/(1-pidach)
      
      # starting values for yi*
      U=rgamma(t,shape=Ni,rate=1+lam0)
      V=rgamma(t,shape=yi,rate=1)
      W=rgamma(t,shape=Ni-yi,rate=lam0)
      yiStar0=-log((U+indi2*W)/(U+indi1*V))
      
      # starting values for regression parameters beta
      if(missing(start)) start=rep(0,dims)
      
      # passing on initial values
      beta_old=beta_new=start
      lam_new=lam0                       
      yiStar_new=yiStar0

      # matrices to save the draws for beta and yi*
      beta=matrix(NA,dims,sim)

      # counter for acceptance rate
      count=numeric(sim)
      
      # calculation of BN for posterior distribution
      trigam=trigamma(Ni)
      Xtilde=X/(2*trigam)
      BN=solve(B0inv+t(Xtilde)%*%X)                             
      BNinv=solve(BN)
      C=chol(BN)
      tC=t(C)
      
      # pre-calculations for bN for posterior distribution
      KN=BN%*%t(Xtilde)
      mN=BN%*%B0invb0
      
      # pre-calculations for acceptance rate
      NiX=Ni*X
      negtwoNi=-2*Ni

      # internal function: calculation of yi* (aggregated utilities)
      yStar=function(lam_new){
                  U=rgamma(t,shape=Ni,rate=1+lam_new)
                  V=rgamma(t,shape=yi,rate=1)
                  W=rgamma(t,shape=Ni-yi,rate=lam_new)
                  -log((U+indi2*W)/(U+indi1*V))
            }
        
      # internal function: calculation of log acceptance rate
      logprob=function(beta_new,beta_old){
                  dbeta=beta_new-beta_old
                  tdbeta=t(dbeta)
                  sumbeta=beta_new+beta_old
                  (sum(NiX%*%dbeta)
                  +sum(negtwoNi*log(1+exp(-yiStar_new+X%*%beta_new)))
                  -sum(negtwoNi*log(1+exp(-yiStar_new+X%*%beta_old)))
                  -1/2*tdbeta%*%B0inv%*%(sumbeta-2*b0)
                  +1/2*tdbeta%*%BNinv%*%(sumbeta-2*bN))            
              }                                                                                         

      # ----------------------- start independence MH-sampler -----------------------
      for(s in seq_len(sim)){
      
          if(s==(burn+1)){stop1=proc.time()[3]}
          
          # Step 1
             
              bN=mN+KN%*%yiStar_new
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
              
          # Step 2    
      		
        	    lam_new=exp(X%*%beta_new)
        	    yiStar_new=yStar(lam_new)
          
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
      out <- list(beta=beta,sim=sim,burn=burn,acc=acc,dims=dims,t=t,b0=b0,B0=B0,duration=duration,
                  duration_wBI=duration_wBI,rate=rate)
      class(out) <- c("binomlogit","binomlogitMH")
      return(out)
}
