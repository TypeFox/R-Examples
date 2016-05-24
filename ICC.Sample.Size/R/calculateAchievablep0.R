calculateAchievablep0<-function(p,k=2,alpha=0.05,tails=2,power=0.80,N)
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
  ##Calulate Fp
  Fp<-(1+(k-1)*p)/(1-p);
  ##Calcuate L0
  L0<-log(Fp)-(zAlpha+qnorm(power))/(sqrt(((k-1)*(N-1)/2*k)));
  ##From L0 calculate p0
  p0<-(exp(L0)-1)/(exp(L0)+k-1);
  ##Create a result dataframe
  resultFrame = data.frame(p0=p0,N=N,p=p,k=k,alpha=alpha,tails=tails,power=power);
  ##Put result dataframe into result list
  result <- list(resultFrame);
  ##Return result list
  return(result)
}