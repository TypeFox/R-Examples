dRUMAuxMix <-
function(yi,Ni,X,sim=12000,burn=2000,b0,B0,start,verbose=500){
  
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
      beta_new=start
      lam_new=lam0                       
      yiStar_new=yiStar0
      
      # matrices to save the draws for beta and yi*
      beta=matrix(NA,dims,sim)
      
      # pre-calculations for BN for posterior distribution
      tX=t(X)
      
      # parameters for the components of the normal mixture distribution
      H=5  # max. number of components
      
      mixture=cbind(matrix(0,nrow=t,ncol=5),matrix(1,nrow=t,ncol=5))
      for(i in seq_len(t)){
        mixture[i,1:length(compmix(Ni[i])$probs)]=compmix(Ni[i])$probs
        mixture[i,6:(6+length(compmix(Ni[i])$probs)-1)]=compmix(Ni[i])$var
      }
      
      probs=mixture[,1:5]             # splitting matrix with probs and vars in two parts
      vars=mixture[,6:10]
      vecvars=as.vector(t(vars))
      st=H*(0:(t-1))
      sd=sqrt(vars)
      logpr=log(probs)
      logsd=log(sd)
      
      # vectors to save chosen component and corresponding variance for each draw
      pos=numeric(t)
      sr2=numeric(t)
      
      # pre-calculations for efficient drawing of scaling factors
      trickmat=outer(seq_len(H),seq_len(H),"<=")         # upper triangular matrix with TRUEs
      vertfkt=matrix(0,nrow=t,ncol=H)                    # matrix for cdfs
      pre=probs/sd
      twovar=-2*vars
      
      # internal function: calculation of yi* (aggregated utilities)
      yStar=function(lam_new){
                  U=rgamma(t,shape=Ni,rate=1+lam_new)
                  V=rgamma(t,shape=yi,rate=1)
                  W=rgamma(t,shape=Ni-yi,rate=lam_new)
                  -log((U+indi2*W)/(U+indi1*V))
            }
      
      # ----------------------- start auxiliary mixture sampler -----------------------
      for(s in seq_len(sim)){
        
        if(s==(burn+1)){stop1=proc.time()[3]}
        
        # Step 1 - Mixture Components
        
            # non-standardized probabilities
            Pr=pre*exp(rep.int((yiStar_new-log(lam_new))^2,H)/twovar)
            # non-standardized cdfs
            vertfkt=Pr%*%trickmat
            # inversion method
            pos=rowSums(vertfkt[,H]*runif(t) > vertfkt)+1
            sr2=vecvars[st+pos]         # current scaling factors
        
        # Step 2 - Betas
        
            Xtilde=X/sr2
            sum1=tX%*%Xtilde
            BN=solve(B0inv+sum1)
            tC=t(chol(BN))
            
            XtYstar=X*yiStar_new
            sum2=colSums(XtYstar/sr2)
            mN=B0invb0+sum2
            bN=BN%*%mN
            beta_new=as.vector(bN+tC%*%rnorm(dims,0,1))
        
        # Step 3 - Utilities
        
            lam_new=exp(X%*%beta_new)
            yiStar_new=yStar(lam_new)
            
        # save new draws
        beta[,s]=beta_new
        
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
      # ----------------------- end auxiliary mixture sampler -----------------------
      
      # stop timer
      finish=proc.time()[3]
      duration=finish-starttime
      duration_wBI=finish-stop1
      
      # output
      out <- list(beta=beta,sim=sim,burn=burn,dims=dims,t=t,b0=b0,B0=B0,duration=duration,duration_wBI=duration_wBI)
      class(out) <- "binomlogit"
      return(out)
}
