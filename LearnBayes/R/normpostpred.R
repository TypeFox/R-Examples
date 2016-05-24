normpostpred=function(parameters,sample.size,f=min)
{
  normalsample=function(j,parameters,sample.size)
    rnorm(sample.size,mean=parameters$mu[j],sd=sqrt(parameters$sigma2[j]))
  
  m=length(parameters$mu)
  post.pred.samples=sapply(1:m,normalsample,parameters,sample.size)

  stat=apply(post.pred.samples,2,f)

  return(stat)
}