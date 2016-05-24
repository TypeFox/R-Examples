#' Performs social network based diffusion analysis in a Bayesian context
#' @param formatteddata  formatted data  
#' @param its number of iterations
#' @param pilot_tuner1 tuning parameter for the social effect
#' @param pilot_tuner2 tuning parameter for the asocial effect
#' @param start1 start value for the social parameter
#' @param start2 start value for the asocial parameter
#' @export

library(mvtnorm)
rjmcmc = function(formatteddata,its,pilot_tuner1,pilot_tuner2,start1,start2,start3,p1,p2){
  
  TimeD = formatteddata[[1]][,7] # time interval of solving
  censored = formatteddata[[1]][,6] # 1/0 binary variable indicating censored status after expt end
  Aij = formatteddata[[1]][,8]# interaction covariate based on network
  NaiveD = formatteddata[[2]] # naive status at each time point for each unique individual
  
  
  
  s0= start1 # strength of social transmission
  baseline_rate = lambda0 = start2# baseline rate of acquisition in absence of social transimission
  model0 = start3
  acceptcounter = 0 #(not used)
  Jumbo<-array(0,c(its,3)) # storage of updated parameters (except interactions) at each iteration
  newparam<-array(0.5,3) # allocation of storage during update process
  CurrentParam<-array(0.5,3)# allocation of storage during update process
  
  newparam[1]=CurrentParam[1]<-s0    # used
  newparam[2]=CurrentParam[2]<-lambda0 # used
  newparam[3]=CurrentParam[3]<-model0 # used
  
  #----------------------------------
  #             FUNCTIONS
  #----------------------------------
  
  
  # updating the s parameter (allowed to take only positive values)
  
  updates<-function(CurrentParam,newparam,model){
    GU3<-runif(1,CurrentParam[1]-pilot_tuner1,CurrentParam[1]+ pilot_tuner1)
    proposal = c(GU3,CurrentParam[2])
    num<-CpS(proposal,model)[[1]]
    den<-CpS(CurrentParam,model)[[1]]
    acc<-exp(num-den)
    acceptr<-min(1,acc)
    r<-runif(1)
    newparam[1]<-ifelse((r<=acceptr),GU3,CurrentParam[1])
    return(newparam[1])
  }
  
  
  # updating the lambda (baseline rate of acquisition) parameter (allowed to take only positive values)
  updatelambda<-function(CurrentParam,newparam,model){
    GU3<-runif(1,CurrentParam[2]-pilot_tuner2,CurrentParam[2]+ pilot_tuner2)      
    proposal = c(CurrentParam[1],GU3)
    num<-CpS(proposal,model)[[1]]
    den<-CpS(CurrentParam,model)[[1]]
    acc<-exp(num-den)
    acceptt<-min(1,acc)
    r<-runif(1)
    newparam[2]<-ifelse((r<=acceptt),GU3,CurrentParam[2])
    acceptcounter<-ifelse((r<=acceptt),1,0)
    list(newparam[2],acceptcounter)
  }
  
  
  
  CpS = function(parameterproposal,model){
    baseline = exp(parameterproposal[2])
    social_rate = exp(parameterproposal[1])
    
    if(model==1){hazard = baseline}
    if(model==2){hazard = baseline + (social_rate)*Aij}
    
    uncensored = 1-censored
    log_likelihood_u = sum(log(hazard*exp(-hazard*TimeD))*uncensored) + sum(-hazard*TimeD*NaiveD)    
    log_likelihood_c =  sum(-hazard*censored)
    log_likelihood = log_likelihood_u + log_likelihood_c     
    lambdaprior<-  log(dunif(parameterproposal[2],-p1,p1)) # s prior
    
    if(model==1){sprior = 0}
    if(model==2){sprior = log(dunif(parameterproposal[1],-p2,p2))}
    
    pzoid<-log_likelihood + lambdaprior  + sprior
    pzoid
  }
  
  
  
  model_space_travel= function(model,t,CurrentParam,Jumbo){
    
    book_flight<-function(model){
      booking<-sample(c(1,2),1,prob=c(rep(0.5,2)))
      while(booking==model){
        booking<-sample(c(1,2),1,prob=c(rep(0.5,2)))
      }
      return(booking)
    }
    
    flight<-book_flight(model)
    
    
    pkg1= (c(rnorm(n=1,mean1,sd1) ))
    pkg2= (c(rmvnorm(n=1,mean2,sigma2) ))
    
    
    
    pkg=function(flight){
      switch(flight,
{#flight=1
  apple=((pkg1[1]<p1)&&(pkg1[1]>-p1))
  apple
},
             
{#flight=2
  bear=((pkg2[2]<p1)&&(pkg2[2]>-p1)&&(pkg2[1]<p2)&&(pkg2[1]>-p2))
  bear
}
      )
    }
    
    green_light=pkg(flight)
    #green_light = TRUE
    if(green_light==TRUE)  # TRUE implies that the simulated values are non zero
      
      
    {
      
      if (model==1){num_proposal_density=log(dnorm(x=c(CurrentParam[2]),mean1,sd1))}
      if (model==2){num_proposal_density=log(dmvnorm(x=c(CurrentParam[c(1:2)]),mean2,sigma2))}                    
      if(flight==1){den_proposal_density=log(dnorm(x=c(pkg1[1]),mean1,sd1))}
      if(flight==2){den_proposal_density=log(dmvnorm(x=c(pkg2[c(1:2)]),mean2,sigma2))}          
      if(flight==1){num_posterior_density=CpS(c(CurrentParam[1],pkg1[1]),1)[[1]]}
      if(flight==2){num_posterior_density=CpS(c(pkg2[1:2]),2)[[1]]}
      if(model==1){den_posterior_density=CpS(c(CurrentParam[1:2]),1)[[1]]}
      if(model==2){den_posterior_density=CpS(c(CurrentParam[1:2]),2)[[1]]}
      num=    num_posterior_density + num_proposal_density        #   (they are both in logs hence the addition)
      den=    den_posterior_density + den_proposal_density
      A<-exp(num - den)
    }
    
    
    
    if(green_light==FALSE) {A = 0 }
    accept<-min(1,A)
    checking<-runif(1)
    if (checking <= accept)
    {
      model<-flight                	
      Jumbo[t,3] = model              		
      
      if(flight==2){
      Jumbo[t,1]=pkg2[1]
      Jumbo[t,2]= pkg2[2]
      }
      if(flight==1){
        Jumbo[t,1]= 0
        Jumbo[t,2] = pkg1[1]
      }
      
    }else
      if(checking>accept) {
        Jumbo[t,1:3]=CurrentParam[1:3]
        
      }
    
    list(model,Jumbo[t,1],Jumbo[t,2])
    
  }
  
 #---------------------------------------------- 
  
  model = 2
  for(t in 1:its){
    CurrentParam[1]=Jumbo[t,1]=updates(CurrentParam,newparam,model)[[1]] # s
    CurrentParam[2]=Jumbo[t,2]=updatelambda(CurrentParam,newparam,model)[[1]] # lambda
  }
  
  m2 = Jumbo
  burnin = its/10
  
  par(mfrow=c(2,2))
  plot(Jumbo[burnin:its,1],type="l",col="blue",ylab="social effect",main="Trace plot for social effect, s' ",lwd=2)
  plot(Jumbo[burnin:its,2],type="l",col="red",ylab="asocial effect",main="Trace plot for asocial effect, lambda0' ",lwd=2)
  plot(density(Jumbo[burnin:its,1],adjust=3),col="darkblue",main="Density plot of social effect, s'",lwd=3)
  acf(Jumbo[burnin:its,1],main="ACF plot for social effect, s'")
  
  
  model = 1
  for(t in 1:its){
    CurrentParam[2]=Jumbo[t,2]=updatelambda(CurrentParam,newparam,model)[[1]] # lambda
  }
  m1 = Jumbo
  
  
  mean1 = mean(m1[,2][burnin:its])
  sd1 = sd(m1[,2][burnin:its])
  mod1 = rnorm(1,mean1,sd1)
  
  mb11=mean(m2[,1][burnin:its])
  sdb11=sd(m2[,1][burnin:its])^2
  mb21=mean(m2[,2][burnin:its])
  sdb21=sd(m2[,2][burnin:its])^2
  
  b1b2=cov(m2[,1][burnin:its],m2[,2][burnin:its]) # the first column has the iteration number
  sigma2= matrix(c(sdb11,b1b2,b1b2,sdb21),ncol=2)
  mean2=c(mb11,mb21)  
  
  
  
 #------------------------------------- 
  
  
  model=model0
  for(t in 1:its){
    switch(model,
           
{     #m1
  CurrentParam[2]=Jumbo[t,2]=newparam[2]=updatelambda(CurrentParam,newparam,model)[[1]]
  
},
           
{    #m2
  
  CurrentParam[1]=Jumbo[t,1]=updates(CurrentParam,newparam,model)[[1]] # s
  CurrentParam[2]=Jumbo[t,2]=updatelambda(CurrentParam,newparam,model)[[1]] # lambda
  
  
}
    )
    travel<-model_space_travel(model,t,CurrentParam,Jumbo)
    model = CurrentParam[3]=Jumbo[t,3]<-travel[[1]]
    CurrentParam[1]=Jumbo[t,1]<-travel[[2]]
    CurrentParam[2]=Jumbo[t,2]<-travel[[3]]
    
  }
  burnin = its/10
  
 #plot(c(burnin:its),(Jumbo[burnin:its,3]),type="l",col="blue",ylab="Model",main="Model trace plot ")
  models = summary(as.factor(Jumbo[burnin:its,3]))    
  rjmcmcresults = list(models)  
  rjmcmcresults
  
  
}