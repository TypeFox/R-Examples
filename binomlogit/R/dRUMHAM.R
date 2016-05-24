dRUMHAM <-
function(yi,Ni,X,sim=12000,burn=2000,b0,B0,start,low=0.05,up=0.95,verbose=500){
  
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
      
      # matrix to save the draws for beta
      beta=matrix(NA,dims,sim)
      
      # counter for acceptance rate
      count=numeric(sim)
            
      # indicators distinguishing between observations used for AuxMix/IndMH and pre-calculations
      indA=as.numeric((yi/Ni)<=low|(yi/Ni)>=up)    # AuxMix
      indApos=which(indA==1)
      indAcount=length(indApos)
      indAsr2=numeric(t)
      XindA=X*indA

      indS=as.numeric((yi/Ni)>low&(yi/Ni)<up)      # MH
      indSpos=which(indS==1)
      indScount=length(indSpos)
      XindS=X*indS
      
      # pre-calculations for BN for posterior distribution
      tX=t(X)
      trigam=2*trigamma(Ni)
      trigamS=indS/trigam
      BSinv=B0inv+tX%*%(X*trigamS)
      BS=solve(BSinv)
      
      # pre-calculations for bN for posterior distribution
      partsum=t(XindS)/trigam
      
      # pre-calculations for acceptance rate
      NiX=Ni*X
      NiXS=Ni*XindS
      negtwoNis=-2*Ni*indS
      negtwoNi=-2*Ni[indSpos]
      trigamSpos=trigam[indSpos]
      
      # parameters for the components of the normal mixture distribution
      H=5  # max. number of components
      
      mixture=cbind(matrix(0,nrow=t,ncol=5),matrix(1,nrow=t,ncol=5))
      for(i in which(indA==1)){
        mixture[i,1:length(compmix(Ni[i])$probs)]=compmix(Ni[i])$probs
        mixture[i,6:(6+length(compmix(Ni[i])$probs)-1)]=compmix(Ni[i])$var
      }
      mixture=matrix(mixture[indApos,],nrow=indAcount)
      
      probs=mixture[,1:5]               # splitting matrix with probs and vars in two parts
      vars=mixture[,6:10]
      vecvars=as.vector(t(vars))
      st=H*(0:(indAcount-1))
      sd=sqrt(vars)
      logpr=log(probs)
      logsd=log(sd)
      
      # vectors to save chosen component and corresponding variance for each draw
      pos=numeric(indAcount)
      sr2=numeric(indAcount)
      
      # pre-calculations for efficient drawing of scaling factors
      trickmat=outer(seq_len(H),seq_len(H),"<=")         # upper triangular matrix with TRUEs
      vertfkt=matrix(0,nrow=indAcount,ncol=H)            # matrix for cdfs
      pre=probs/sd
      twovar=-2*vars
      
      # internal function: calculation of yi* (aggregated utilities)
      yStar=function(lam_new){
                  U=rgamma(t,shape=Ni,rate=1+lam_new)
                  V=rgamma(t,shape=yi,rate=1)
                  W=rgamma(t,shape=Ni-yi,rate=lam_new)
                  -log((U+indi2*W)/(U+indi1*V))
            }
      
      # internal function: calculation of log acceptance rate
      logprob=function(beta_new,beta_old,XtYstar){
                  mS=B0invb0+colSums(XtYstar*trigamS)
                  bS=BS%*%mS
                  dbeta=beta_new-beta_old
                  tdbeta=t(dbeta)
                  sumbeta=beta_new+beta_old
                  (sum(NiXS%*%dbeta)
                  +sum(negtwoNis*log(1+exp(-yiStar_new+X%*%beta_new)))
                  -sum(negtwoNis*log(1+exp(-yiStar_new+X%*%beta_old)))
                  -1/2*tdbeta%*%B0inv%*%(sumbeta-2*b0)
                  +1/2*tdbeta%*%BSinv%*%(sumbeta-2*bS))         
            }
      
      # -------------------------- start HAM sampler --------------------------
      for(s in seq_len(sim)){
        
          if(s==(burn+1)){stop1=proc.time()[3]}
          
          # Step 1 - Scaling Factors für AuxMix
          
              # non-standardized probabilities
              Pr=pre*exp(rep.int((yiStar_new[indApos]-log(lam_new[indApos]))^2,H)/twovar)
              # non-standardized cdfs
              vertfkt=Pr%*%trickmat
              # inversion method
              pos=rowSums(vertfkt[,H]*runif(indAcount) > vertfkt)+1
              sr2=vecvars[st+pos]                  # current scaling factors
              indAsr2[indApos]=1/sr2
          
          # Step 2 - Betas
          
              sum1=matrix(tX[,indApos],nrow=dims)%*%matrix(X[indApos,]/sr2,ncol=dims)
              BN=solve(BSinv+sum1)
              tC=t(chol(BN))
              
              XtYstar=X*yiStar_new
              sum2=colSums(XtYstar*(trigamS+indAsr2))
              mN=B0invb0+sum2               
              bN=BN%*%mN
              
              beta_new=as.vector(bN+tC%*%rnorm(dims,0,1))
              if(indScount>0){
              u=runif(1,0,1)
              prob=exp(logprob(beta_new,beta_old,XtYstar))                           
              if(u<prob){
              count[s]=1
              } else {
              count[s]=0
              beta_new=beta_old
              }}
          
          # Step 3 - Utilities
          
              lam_new=exp(X%*%beta_new)
              yiStar_new=yStar(lam_new)
        
          # Abspeichern der neuen Samples
          beta[,s]=beta_new
          beta_old=beta_new
          
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
      # -------------------------- end HAM sampler --------------------------
      
      # stop timer
      finish=proc.time()[3]
      duration=finish-starttime
      duration_wBI=finish-stop1
      
      # acceptance rate
      rate=round(sum(count[(burn+1):sim])/(sim-burn)*100,2)
      
      # output
      out <- list(beta=beta,sim=sim,burn=burn,dims=dims,t=t,b0=b0,B0=B0,low=low,up=up,duration=duration,
                  duration_wBI=duration_wBI,rate=rate)
      class(out) <- c("binomlogit","binomlogitHAM")
      return(out)
}
