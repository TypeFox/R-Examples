calculateIccSampleSize<-function(p=0,p0=0,k=2,alpha=0.05,tails=2,power=0.80,by="",step=0.05)
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
  ##Calculate Fp and Fp0
  Fp<-(1+(k-1)*p)/(1-p);
  Fp0<-(1+(k-1)*p0)/(1-p0)
  ##Calculate N, rounded up to nearest integer
  Nraw<-(1+(2*(zAlpha+qnorm(power))^2*k)/((log(Fp/Fp0))^2*(k-1)));
  N<-ceiling(Nraw)
  ##Place results in a datafrome
  resultdf <- data.frame(N=N,p=p,p0=p0,k=k,alpha=alpha,tails=tails,power=power)
  ##If by is empty then result is resultdf
  if(by == "")
  {
    result<-list(resultdf)
  }
  ##If by is "p", then calculate sample size for every p at intervals of step between 0 and 1
  if(by == "p")
  {
    ##Create blank vectors to store p's and N's
    pVector<-c();
    nVector<-c();
    ##Zero p
    p<-0;
    ##Calculate sample size at each step
    while(p<=1)
    {
      ##Calculate N
      Fp<-(1+(k-1)*p)/(1-p);
      Fp0<-(1+(k-1)*p0)/(1-p0);
      Nraw<-(1+(2*(zAlpha+qnorm(power))^2*k)/((log(Fp/Fp0))^2*(k-1)));
      N<-ceiling(Nraw);
      ##Add N and p to vector of N and p
      nVector<-c(nVector,N);
      pVector<-c(pVector,p);
      ##Set p for next step of loop
      p<-p+step
    }
    ##Create dataframe from lists of p and N
    sampleSize <- data.frame(p=pVector,N=nVector);
    ##Make results from resultdf and samplesize
    result<-list(resultdf,sampleSize)
  }
  ##If by is "p0", then calculate sample size for every p at intervals of step between 0 and 1
  if(by == "p0")
  {
    ##Create blank vectors to store p0's and N's
    p0Vector<-c();
    nVector<-c();
    ##Zero p0
    p0<-0;
    ##Calculate sample size at each step
    while(p0<=1)
    {
      ##Calculate N
      Fp<-(1+(k-1)*p)/(1-p);
      Fp0<-(1+(k-1)*p0)/(1-p0);
      Nraw<-(1+(2*(zAlpha+qnorm(power))^2*k)/((log(Fp/Fp0))^2*(k-1)));
      N<-ceiling(Nraw);
      ##Add N and p0 to vector of N and p0
      nVector<-c(nVector,N);
      p0Vector<-c(p0Vector,p0);
      ##Set p0 for next step of loop
      p0<-p0+step
    }
    ##Create dataframe from lists of p and N
    sampleSize <- data.frame(p0=p0Vector,N=nVector);
    ##Make results from resultdf and samplesize
    result <- list(resultdf,sampleSize)
  }
  ##If by is "both", then calculate sample size for every combination of p and p0 at intervals of step between 0 and 1 
  if(by == "both")
  {
    ##Create blank vectors to store N's
    nVector<-c();
    ##Zero p and p0
    p<-0;
    p0<-0;
    ##Create a counter to measure repeats of p0 loop
    p0Count<-0;
    while(p0<=1)
    {
      ##Create a counter to measure repeats of p loop
      pCount<-0;
      ##Calculate sample size at each step of p
      while(p<=1)
      {
        ##Calculate N
        Fp<-(1+(k-1)*p)/(1-p);
        Fp0<-(1+(k-1)*p0)/(1-p0);
        Nraw<-(1+(2*(zAlpha+qnorm(power))^2*k)/((log(Fp/Fp0))^2*(k-1)));
        N<-ceiling(Nraw);
        ##Add N to vector of N's
        nVector<-c(nVector,N);
        ##Set p for next step of p loop
        p<-p+step;
        ##Add 1 to counter of p loop repeats
        pCount<-pCount+1
      }
      ##Zero p for next step of p0 loop
      p<-0;
      ##Set p0 for next step of p0 loop
      p0<-p0+step;
      ##Add 1 to counter of p0 loop repeats 
      p0Count<-p0Count+1
    }
    ##Create p axis labels
    ##Create a counter
    pLabelCount<-0;
    ##Create a vector with p axis labels
    pLabelVector<-c();
    ##Add labels for each each p
    while(pLabelCount<pCount)
    {
      pLabelVector<-c(pLabelVector,step*pLabelCount);
      pLabelCount<-pLabelCount+1
    }
    ##Create p0 axis labels
    ##Create a counter
    p0LabelCount<-0;
    ##Create a vector with p0 axis labels
    p0LabelVector<-c();
    ##Add labels for each p0
    while(p0LabelCount<p0Count)
    {
      p0LabelVector<-c(p0LabelVector,step*p0LabelCount);
      p0LabelCount<-p0LabelCount+1
    }
    ##Create a matrix from vector of N's
    nMatrix <- matrix(nVector,nrow=pCount,ncol=p0Count);
    ##Create a dataframe from the above matrix
    nDataframe<-as.data.frame(nMatrix);
    ##Add labels for p and p0
    row.names(nDataframe)<-pLabelVector;
    colnames(nDataframe)<-p0LabelVector;
    ##Make results from resultdf and the dataframe above
    result<-list(resultdf,nDataframe)
  }
  ##return the results
  return(result)
}
