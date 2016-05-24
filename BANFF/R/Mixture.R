Mixture=function(data,pre)
{
  resultall=0
  for(i in 1:length(pre$mean))
  {
    resultall=resultall+pre$pro[i]/sum(pre$pro)*dnorm(data,mean=pre$mean[i],sd=pre$sd[i])
  }
  return(resultall)
}
