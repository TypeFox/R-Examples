#' Performs social network based diffusion analysis in a Bayesian context
#' @param formatteddata  formatted data  
#' @param its number of iterations
#' @param pilot_tuner1 tuning parameter for the social effect
#' @param pilot_tuner2 tuning parameter for the asocial effect
#' @param start1 start value for the social parameter
#' @param start2 start value for the asocial parameter
#' @export


mcmcre = function(formatteddata,its,pilot_tuner1,pilot_tuner2,start1,start2){
  
  TimeD = formatteddata[[1]][,7] # time interval of solving
  censored = formatteddata[[1]][,6] # 1/0 binary variable indicating censored status after expt end
  Aij = formatteddata[[1]][,8]# interaction covariate based on network
  NaiveD = formatteddata[[2]] # naive status at each time point for each unique individual
  aptitudes = formatteddata[[1]][,11] # random effects at the individual level
  laptitudes = length(aptitudes)
  tuner = c(pilot_tuner1,pilot_tuner2,0.01,rep(0.05,laptitudes))
  
  
  s0= start1 # strength of social transmission
  baseline_rate = lambda0 = start2# baseline rate of acquisition in absence of social transimission
  gamma = array(0,c(length(aptitudes),1))
  numparam = 3 + length(aptitudes)# social, baseline, variances, and aptitudes
  acceptcounter = 0 #(not used)
  Jumbo<-array(0,c(its,numparam)) # storage of updated parameters (except interactions) at each iteration
  newparam<-array(0.5,numparam) # allocation of storage during update process
  CurrentParam<-array(0.5,numparam)# allocation of storage during update process
  
  newparam[1]=CurrentParam[1]<-s0    # used
  newparam[2]=CurrentParam[2]<-lambda0 # used
  newparam[3:numparam]=CurrentParam[3:numparam]<-0.1 # used
  
  
  #----------------------------------
  #             FUNCTIONS
  #----------------------------------
  
  
  
  blockupdate<-function(CurrentParam){
    block = CurrentParam
    for(i in 1:numparam){
      block[i] = runif(1,CurrentParam[i]-tuner[i],CurrentParam[i]+tuner[i])
    }
    newparam = CurrentParam 
    num<-CpS(block)
    den<-CpS(CurrentParam)
    acc<-exp(num-den)   
    acceptr<-min(1,acc)
    r<-runif(1)
    
    if(r<=acceptr){newparam = block}
    if(r>acceptr){newparam = CurrentParam}       
    return(newparam)
    
  }
  
   
  
  
  CpS = function(parameterproposal){
    baseline = exp(parameterproposal[2])
    social_rate = exp(parameterproposal[1])
    smartness = c(exp(parameterproposal[4:numparam]))
    #smartness = 1
    hazard = (baseline * smartness) + (social_rate)*Aij # hazard function
    uncensored = 1-censored
    log_likelihood_u = sum(log(hazard*exp(-hazard*TimeD))*uncensored) + sum(-hazard*TimeD*NaiveD)    
    log_likelihood_c =  sum(-hazard*censored)
    log_likelihood = log_likelihood_u + log_likelihood_c     
    lambdaprior<-  log(dunif(parameterproposal[2],-10,10)) # s prior
    sprior<- log(dunif(parameterproposal[1],-10,10))# lambda prior
    varianceprior = log(dnorm(parameterproposal[3],0.1))
    rsd = sqrt(exp(parameterproposal[3]))
    smartnessprior<- sum(log(dnorm(parameterproposal[4:numparam],0,rsd)))# lambda prior
    pzoid<-log_likelihood + lambdaprior  + sprior + smartnessprior + varianceprior
    pzoid
  }
  
  
  
  for(t in 1:its){
   
    CurrentParam = Jumbo[t,] =blockupdate(CurrentParam)
   
    }
  
  burnin = its/10
  
  
  par(mfrow=c(2,2))
  plot(Jumbo[burnin:its,1],type="l",col="blue",ylab="social effect",main="Trace plot for social effect, s' ",lwd=2)
  plot(Jumbo[burnin:its,2],type="l",col="red",ylab="asocial effect",main="Trace plot for asocial effect, lambda0' ",lwd=2)
  plot(density(Jumbo[burnin:its,1],adjust=3),col="darkblue",main="Density plot of social effect, s'",lwd=3)
  acf(Jumbo[burnin:its,1],main="ACF plot for social effect, s'")
  
  
  reffectsrows = 3*(numparam - 3)
  restore=matrix(0,nrow = reffectsrows,ncol = 2)
  my_re_table=matrix(0,nrow = reffectsrows,ncol = 2)
  colnames(my_re_table)=c("summary","random effects (aptitudes)")
  
  
  for(i in 1:laptitudes){
    start = ifelse(i!=1,1+(3*(i-1)),1)
    end = (3*i)
    columnindex = i
    meanval = mean(Jumbo[,columnindex])
    ci1val =  quantile(Jumbo[,columnindex],0.025)
    ci2val =  quantile(Jumbo[,columnindex],0.975)
    my_re_table[start:end,2]= c(meanval,ci1val,ci2val)}
  
    restore=(my_re_table)
    restore[,2]=round(my_re_table[,2],digits=5)
    restore[1:reffectsrows,1]=c("mean","95%ci1","95%ci2")
    
  for(i in 1:reffectsrows){
   for(j in 1:2){
      if (restore[i,j]==0) {restore[i,j]=""}
    }
  } 
  
 #----------------------- 
  datahouse=matrix(0,9,2)
  my_summary_table=matrix(0,9,2)
  colnames(my_summary_table)=c("summary","model parameters")
  rownames(my_summary_table)=c("lambda0","","","s","","","sigma2","","")
  
  #--------------model 1 Summary--------------------------------------------------
  
  # lambda
  my_summary_table[1:3,2]=
    c(mean(Jumbo[burnin:its,2]),
      quantile(Jumbo[burnin:its,2],0.025),
      quantile(Jumbo[burnin:its,2],0.975)
      
    )
  
  #------------- model 2 Summary--------------------------------------------------
  
  # s
  my_summary_table[4:6,2]=
    c(mean(Jumbo[burnin:its,1]),
      quantile(Jumbo[burnin:its,1],0.025),
      quantile(Jumbo[burnin:its,1],0.975)
      
    )
  
  #sigma2
  my_summary_table[7:9,2]=
    c(mean(Jumbo[burnin:its,3]),
      quantile(Jumbo[burnin:its,3],0.025),
      quantile(Jumbo[burnin:its,3],0.975)
      
    )
  
  
  datahouse=(my_summary_table)
  datahouse[,2]=round(my_summary_table[,2],digits=5)
  datahouse[1:9,1]=
    c("mean",
      "95%ci1",
      "95%ci2"
    )
  #
  
  # ----------------changing all the zeros to ""'s -------------------------------
  
  for(i in 1:9){
    for(j in 1:2){
      if (datahouse[i,j]==0) {datahouse[i,j]=""}
    }
  }
  
  
  list(Jumbo,restore,datahouse)
}
  
  
  
  
