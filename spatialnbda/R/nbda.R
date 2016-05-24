
#' Performs social network based diffusion analysis in a Bayesian context
#' @param formatteddata  formatted data  
#' @param its number of iterations
#' @param pilot_tuner1 tuning parameter for the social effect
#' @param pilot_tuner2 tuning parameter for the asocial effect
#' @param start1 start value for the social parameter
#' @param start2 start value for the asocial parameter
#' @export


  mcmc = function(formatteddata,its,pilot_tuner1,pilot_tuner2,start1,start2){
    
    TimeD = formatteddata[[1]][,7] # time interval of solving
    censored = formatteddata[[1]][,6] # 1/0 binary variable indicating censored status after expt end
    Aij = formatteddata[[1]][,8]# interaction covariate based on network
    NaiveD = formatteddata[[2]] # naive status at each time point for each unique individual
    
    
    
    
    
    s0= start1 # strength of social transmission
    baseline_rate = lambda0 = start2# baseline rate of acquisition in absence of social transimission
    acceptcounter = 0 #(not used)
    Jumbo<-array(0,c(its,3)) # storage of updated parameters (except interactions) at each iteration
    newparam<-array(0.5,2) # allocation of storage during update process
    CurrentParam<-array(0.5,2)# allocation of storage during update process
    
    newparam[1]=CurrentParam[1]<-s0    # used
    newparam[2]=CurrentParam[2]<-lambda0 # used
    
    
    #----------------------------------
    #             FUNCTIONS
    #----------------------------------
    
    # updating the s parameter (allowed to take only positive values)
    
    blockupdate<-function(CurrentParam,newparam){
      GU1<-runif(1,CurrentParam[1]-pilot_tuner1,CurrentParam[1]+ pilot_tuner1)
      GU2<-runif(1,CurrentParam[2]-pilot_tuner2,CurrentParam[2]+ pilot_tuner2)
      proposal = c(GU1,GU2)
      num<-CpS(proposal)[[1]]
      den<-CpS(CurrentParam)[[1]]
      acc<-exp(num-den)
      acceptr<-min(1,acc)
      r<-runif(1)
      newparam[1]<-ifelse((r<=acceptr),GU1,CurrentParam[1])
      newparam[2]<-ifelse((r<=acceptr),GU2,CurrentParam[2])
      list(newparam[1],newparam[2])
    }
    
   
    
    CpS = function(parameterproposal){
      baseline = exp(parameterproposal[2])
      social_rate = exp(parameterproposal[1])
      hazard = baseline + (social_rate)*Aij # hazard function
      uncensored = 1-censored
      log_likelihood_u = sum(log(hazard*exp(-hazard*TimeD))*uncensored) + sum(-hazard*TimeD*NaiveD)    
      log_likelihood_c =  sum(-hazard*censored)
      log_likelihood = log_likelihood_u + log_likelihood_c     
      lambdaprior<-  log(dunif(parameterproposal[2],-10,10)) # s prior
      sprior<- log(dunif(parameterproposal[1],-10,10))# lambda prior
      pzoid<-log_likelihood + lambdaprior  + sprior
      dhatcomponent = -2*log_likelihood
      list(pzoid,(dhatcomponent))
    }
    
    
   
    for(t in 1:its){
      runmodel = blockupdate(CurrentParam,newparam)
      CurrentParam[2]=Jumbo[t,2]= runmodel[[2]] # lambda
      CurrentParam[1]=Jumbo[t,1]= runmodel[[1]]
      Jumbo[t,3] = CpS(CurrentParam)[[2]]
          
    }
    
    burnin = its/10
    Dhat = mean(Jumbo[burnin:its,3]) # the mean of the deviance
    posteriorparameters = c(mean(Jumbo[burnin:its,1]),mean(Jumbo[burnin:its,2]))
    Phat = CpS(posteriorparameters)[[2]]
    Twopd = 2*(Dhat - Phat)
    DIC = Dhat + Twopd 
    
    par(mfrow=c(2,2))
    plot(Jumbo[burnin:its,1],type="l",col="blue",ylab="social effect",main="Trace plot for social effect, s' ",lwd=2)
    plot(Jumbo[burnin:its,2],type="l",col="red",ylab="asocial effect",main="Trace plot for asocial effect, lambda0' ",lwd=2)
    plot(density(Jumbo[burnin:its,1],adjust=3),col="darkblue",main="Density plot of social effect, s'",lwd=3)
    acf(Jumbo[burnin:its,1],main="ACF plot for social effect, s'")
    
    
    
    
    #params = c(mean(Jumbo[burnin:its,1]),mean(Jumbo[burnin:its,2]))
    #creds = c(sd(Jumbo[burnin:its,1]),sd(Jumbo[burnin:its,2]))
    #mcmcresults = list(Jumbo,params,creds)
    
    #mcmcresults
    
    datahouse=matrix(0,6,2)
    my_summary_table=matrix(0,6,2)
    colnames(my_summary_table)=c("summary","full model")
    rownames(my_summary_table)=c("lambda0","","","s'","","")
    
    # lambda
    my_summary_table[1:3,2]=
      c(mean(Jumbo[burnin:its,2]),
        quantile(Jumbo[burnin:its,2],0.025),
        quantile(Jumbo[burnin:its,2],0.975)
        
      )
    
    #s
    my_summary_table[4:6,2]=
      c(mean(Jumbo[burnin:its,1]),
        quantile(Jumbo[burnin:its,1],0.025),
        quantile(Jumbo[burnin:its,1],0.975)
        
      )
    
    datahouse=(my_summary_table)
    datahouse[,2]=round(my_summary_table[,2],digits=5)
    datahouse[1:6,1]=
      c("mean",
        "95%ci1",
        "95%ci2"
      )
    #
    
    # ----------------changing all the zeros to ""'s -------------------------------
    
    for(i in 1:6){
      for(j in 1:2){
        if (datahouse[i,j]==0) {datahouse[i,j]=""}
      }
    }
    mcmcresults = list(Jumbo,datahouse,DIC)
    mcmcresults
   
    
  }