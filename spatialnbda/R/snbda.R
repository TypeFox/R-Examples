#' Performs social network based diffusion analysis in a Bayesian context
#' @param formatteddata  formatted data  
#' @param its number of iterations
#' @param pilot_tuner1 tuning parameter for the social effect
#' @param pilot_tuner2 tuning parameter for the asocial effect
#' @param start1 start value for the social parameter
#' @param start2 start value for the asocial parameter
#' @export


smcmc = function(formatteddata,its,pilot_tuner1,pilot_tuner2,pilot_tuner3,start1,start2,start3){
  
  TimeD = formatteddata[[1]][,7] # time interval of solving
  censored = formatteddata[[1]][,6] # 1/0 binary variable indicating censored status after expt end
  Aij = formatteddata[[1]][,8]# interaction covariate based on network
  NaiveD = formatteddata[[2]] # naive status at each time point for each unique individual
  spatcov = formatteddata[[1]][,10]# spatial/environmental covariate
  
   
  
  s0= start1 # strength of social transmission
  baseline_rate = lambda0 = start2# baseline rate of acquisition in absence of social transimission
  environmental = beta0 = start3 # spatial covariate
  
  acceptcounter = 0 #(not used)
  Jumbo<-array(0,c(its,3)) # storage of updated parameters (except interactions) at each iteration
  newparam<-array(0.5,3) # allocation of storage during update process
  CurrentParam<-array(0.5,3)# allocation of storage during update process
  
  newparam[1]=CurrentParam[1]<-s0    # used
  newparam[2]=CurrentParam[2]<-lambda0 # used
  newparam[3]=CurrentParam[3]<-beta0 # used
  
  #----------------------------------
  #             FUNCTIONS
  #----------------------------------
  
  
  # updating the s parameter (allowed to take only positive values)
  
  updates<-function(CurrentParam,newparam){
    GU3<-runif(1,CurrentParam[1]-pilot_tuner1,CurrentParam[1]+ pilot_tuner1)
    proposal = c(GU3,CurrentParam[2:3])
    num<-CpS(proposal)[[1]]
    den<-CpS(CurrentParam)[[1]]
    acc<-exp(num-den)
    acceptr<-min(1,acc)
    r<-runif(1)
    newparam[1]<-ifelse((r<=acceptr),GU3,CurrentParam[1])
    return(newparam[1])
  }
  
  
  # updating the lambda (baseline rate of acquisition) parameter (allowed to take only positive values)
  updatelambda<-function(CurrentParam,newparam){
    GU3<-runif(1,CurrentParam[2]-pilot_tuner2,CurrentParam[2]+ pilot_tuner2)      
    proposal = c(CurrentParam[1],GU3,CurrentParam[3])
    num<-CpS(proposal)[[1]]
    den<-CpS(CurrentParam)[[1]]
    acc<-exp(num-den)
    acceptt<-min(1,acc)
    r<-runif(1)
    newparam[2]<-ifelse((r<=acceptt),GU3,CurrentParam[2])
    acceptcounter<-ifelse((r<=acceptt),1,0)
    list(newparam[2],acceptcounter)
  }
  
  
  updatecovariate<-function(CurrentParam,newparam){
    GU3<-runif(1,CurrentParam[3]-5,CurrentParam[3]+ 5)
    proposal = c(CurrentParam[1:2],GU3)
    num<-CpS(proposal)[[1]]
    den<-CpS(CurrentParam)[[1]]
    acc<-exp(num-den)
    acceptr<-min(1,acc)
    r<-runif(1)
    newparam[3]<-ifelse((r<=acceptr),GU3,CurrentParam[3])
    return(newparam[3])
  }
  
    
  
  CpS = function(parameterproposal){
    baseline = exp(parameterproposal[2])
    social_rate = exp(parameterproposal[1])
    spatialC = exp(parameterproposal[3])
    hazard = baseline*exp(spatialC) + (social_rate)*Aij # hazard function
    uncensored = 1-censored
    log_likelihood_u = sum(log(hazard*exp(-hazard*TimeD))*uncensored) + sum(-hazard*TimeD*NaiveD)    
    log_likelihood_c =  sum(-hazard*censored)
    log_likelihood = log_likelihood_u + log_likelihood_c     
    lambdaprior<-  log(dunif(parameterproposal[2],-10,10)) # s prior
    sprior<- log(dunif(parameterproposal[1],-10,10))# lambda prior
    sCprior<- log(dunif(parameterproposal[3],-10,10))# lambda prior
    pzoid<-log_likelihood + lambdaprior  + sprior + sCprior
    pzoid
  }
  
  
  
  for(t in 1:its){
    CurrentParam[1]=Jumbo[t,1]=updates(CurrentParam,newparam)[[1]] # s
    CurrentParam[2]=Jumbo[t,2]=updatelambda(CurrentParam,newparam)[[1]] # lambda
    CurrentParam[3]=Jumbo[t,3]=updatecovariate(CurrentParam,newparam)[[1]] # beta
  }
  
  burnin = its/10
  
  
  par(mfrow=c(2,2))
  plot(Jumbo[burnin:its,1],type="l",col="blue",ylab="social effect",main="Trace plot for social effect, s' ",lwd=2)
  plot(Jumbo[burnin:its,2],type="l",col="red",ylab="asocial effect",main="Trace plot for asocial effect, lambda0' ",lwd=2)
  plot(Jumbo[burnin:its,3],type="l",col="lightgoldenrod",ylab="spatial effect",
       main="Trace plot for spatial effect, beta0' ",lwd=2)
  
  
  
  
  params = c(mean(Jumbo[burnin:its,1]),mean(Jumbo[burnin:its,2]),mean(Jumbo[burnin:its,3]))
  creds = c(sd(Jumbo[burnin:its,1]),sd(Jumbo[burnin:its,2]),sd(Jumbo[burnin:its,3]))
  mcmcresults = list(Jumbo,params,creds)
  
  mcmcresults
  
  
}