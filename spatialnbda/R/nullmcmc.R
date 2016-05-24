#' Performs social network based diffusion analysis in a Bayesian context
#' @param formatteddata  formatted data  
#' @param its number of iterations
#' @param pilot_tuner tuning parameter for the asocial effect
#' @param start start value for the asocial parameter
#' @export


nullmcmc = function(formatteddata,its,pilot_tuner,start){
  
  TimeD = formatteddata[[1]][,7] # time interval of solving
  censored = formatteddata[[1]][,6] # 1/0 binary variable indicating censored status after expt end
  Aij = formatteddata[[1]][,8]# interaction covariate based on network
  NaiveD = formatteddata[[2]] # naive status at each time point for each unique individual
  
  
  
  
  
  
  baseline_rate = lambda0 = start# baseline rate of acquisition in absence of social transimission
  acceptcounter = 0 #(not used)
  Jumbo<-array(0,c(its,3)) # storage of updated parameters (except interactions) at each iteration
  newparam<-array(0.5,2) # allocation of storage during update process
  CurrentParam<-array(0.5,2)# allocation of storage during update process
  
  newparam[1]=CurrentParam[1]<-0    # used
  newparam[2]=CurrentParam[2]<-lambda0 # used
  
  
  #----------------------------------
  #             FUNCTIONS
  #----------------------------------
  
  
  
  
  
  # updating the lambda (baseline rate of acquisition) parameter (allowed to take only positive values)
  updatelambda<-function(CurrentParam,newparam){
    GU3<-runif(1,CurrentParam[2]-pilot_tuner,CurrentParam[2]+ pilot_tuner)      
    proposal = c(CurrentParam[1],GU3)
    num<-CpS(proposal)[[1]]
    den<-CpS(CurrentParam)[[1]]
    acc<-exp(num-den)
    acceptt<-min(1,acc)
    r<-runif(1)
    newparam[2]<-ifelse((r<=acceptt),GU3,CurrentParam[2])
    acceptcounter<-ifelse((r<=acceptt),1,0)
    list(newparam[2],acceptcounter)
  }
  
  
  
  
  
  CpS = function(parameterproposal){
    baseline = exp(parameterproposal[2])
    hazard = baseline
    uncensored = 1-censored
    log_likelihood_u = sum(log(hazard*exp(-hazard*TimeD))*uncensored) + sum(-hazard*TimeD*NaiveD)    
    log_likelihood_c =  sum(-hazard*censored)
    log_likelihood = log_likelihood_u + log_likelihood_c     
    lambdaprior<-  log(dunif(parameterproposal[2],-10,10)) # s prior
    sprior<- 0
    pzoid<-log_likelihood + lambdaprior  + sprior
    dhatcomponent = -2*log_likelihood
    list(pzoid,(dhatcomponent))
  }
  
  
  
  for(t in 1:its){
    
    runmodel = updatelambda(CurrentParam,newparam)
    CurrentParam[2]=Jumbo[t,2]= runmodel[[1]] # lambda
    Jumbo[t,3] = CpS(CurrentParam)[[2]]
  }
  
  burnin = its/10
  
  
  Dhat = mean(Jumbo[burnin:its,3]) # the mean of the deviance
  posteriorparameters = c(mean(Jumbo[burnin:its,1]),mean(Jumbo[burnin:its,2]))
  Phat = CpS(posteriorparameters)[[2]]
  Twopd = 2*(Dhat - Phat)
  DIC = Dhat + Twopd 
   
  plot(Jumbo[burnin:its,2],type="l",col="red",ylab="asocial effect",main="Trace plot for asocial effect, lambda0' ",lwd=2)
  plot(density(Jumbo[burnin:its,2],adjust=3),col="darkblue",main="Density plot of asocial effect, lambda0'",lwd=3)
  acf(Jumbo[burnin:its,2],main="ACF plot for baseline rate, lambda0")
  
  
  
  
  #params = mean(Jumbo[burnin:its,2])
  #creds = c(sd(Jumbo[burnin:its,2]))
  #mcmcresults = list(Jumbo,params,creds)
  
  #mcmcresults
  
  datahouse=matrix(0,3,2)
  my_summary_table=matrix(0,3,2)
  colnames(my_summary_table)=c("summary","null model")
  rownames(my_summary_table)=c("lambda0","","")
  #--------------model 1 Summary--------------------------------------------------
  # lambda
  my_summary_table[1:3,2]=
    c(mean(Jumbo[burnin:its,2]),
      quantile(Jumbo[burnin:its,2],0.025),
      quantile(Jumbo[burnin:its,2],0.975)
    )
  datahouse=(my_summary_table)
  datahouse[,2]=round(my_summary_table[,2],digits=5)
  datahouse[1:3,1]=
    c("mean",
      "95%ci1",
      "95%ci2"
    )
  # ----------------changing all the zeros to ""'s -------------------------------
  for(i in 1:3){
    for(j in 1:2){
      if (datahouse[i,j]==0) {datahouse[i,j]=""}
    }
  }
  
  mcmcresults = list(Jumbo,datahouse,DIC)
  mcmcresults
  
  
}