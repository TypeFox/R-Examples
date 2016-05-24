calculateIccPower<-function(p,p0=0,k=2,alpha=0.05,tails=2,N,by="",desiredPower=0.80,maxN=N*10,step)
{
  ## Assign Z score of alpha depending on number of tails
  if(tails==2)
  {
    zAlpha<-qnorm(1-alpha/2)
  }
  if(tails==1)
  {
    zAlpha<-qnorm(1-alpha)
  }
  ##Caculate Fp and Fp0
  Fp<-(1+(k-1)*p)/(1-p)
  Fp0<-(1+(k-1)*p0)/(1-p0)
  ##Calculate zB
  zB<-sqrt(((k-1)*(N-1))/(2*k))*log(Fp/Fp0)-zAlpha
  ##Calculate power
  power<-pnorm(zB)
  ##Assign power and N vectors
  powerVector<-c(power)
  NVector<-c(N)
  ##If by=N increase by N until desired power or max N is reached
  if(by=="N")
  {
    while(power<=desiredPower&N<=maxN)
    {
      ##If step is missing give it value of 1
      if(missing(step))
        step<-1
      ##Increase N by step
      N<-N+step
      ##Add new N to vector of N's
      NVector<-c(NVector,N)
      ##Caculate Fp and Fp0
      Fp<-(1+(k-1)*p)/(1-p)
      Fp0<-(1+(k-1)*p0)/(1-p0)
      ##Calculate zB
      zB<-sqrt(((k-1)*(N-1))/(2*k))*log(Fp/Fp0)-zAlpha
      ##Calculate power
      power<-pnorm(zB)
      ##Add new power to the vector of powers
      powerVector<-c(powerVector,power)
    }
  }
  ##If by=power perform sample size calculation for each power at intervals of the value assigned to step
  ##between power of the study and desired power until desired power or max N is reached
  if(by=="power")
  {
    while(power<=desiredPower&N<=maxN)
    {
      ##If step is missing assign it a value of 0.05
      if(missing(step))
        step<-0.05
      ##Increase power by step
      power<-power+step
      ##Add new power to vector of powers
      powerVector<-c(powerVector,power)
      ##Calculate Fp and Fp0
      Fp<-(1+(k-1)*p)/(1-p)
      Fp0<-(1+(k-1)*p0)/(1-p0)
      ##Calculate sample size rounded up to nearest integer
      Nraw<-(1+(2*(zAlpha+qnorm(power))^2*k)/((log(Fp/Fp0))^2*(k-1)))
      N<-ceiling(Nraw)
      ##Add N to vector of N's
      NVector<-c(NVector,N)
    }
  }
  ##Create Dataframe of N's and Powers
  NPower <- data.frame(N=NVector, Power=powerVector)
  ##Create a Dataframe with the empirical power of the study and the parameters of the study
  parameters <- data.frame(p=p,p0=p0,k=k,alpha=alpha,tails=tails,N=NVector[1],power=powerVector[1])
  ##If additional powers were calculated results is a list with both the dataframe of N's and Powers
  ##and the dataframe of empirical power and parameters. Else results is a list only of the second dataframe
  if(by=="N"|by=="power")
  {
    result <-list(parameters,NPower) 
  } else {
    result <-list(parameters)
  }
  return(result)
}