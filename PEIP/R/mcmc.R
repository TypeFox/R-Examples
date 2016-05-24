mcmc <-
function(alogprior,aloglikelihood,agenerate,alogproposal,m0,niter)
{  
###   Figure out some size information.
### 
n=length(m0)
### 
###  Allocate space for the results.
### 
mout=matrix(0,n,niter)
### 
###  Initialize the starting point.
### 
mout[,1]=m0
current=m0
lMAP=-Inf
mMAP=current
nacc=0
### 
###  The main loop.
###
generate = get(agenerate)
logproposal =  get(alogproposal)
loglikelihood = get(aloglikelihood)
logprior =get(alogprior)

for(k in 2:niter )
  {
###  Generate a candidate model from the previous model.
### 
  candidate=generate(current)
### 
###  Evalate the logarithm of the acceptance ratio.
### 
  lpcandidate=logprior(candidate)
  
  llcandidate=loglikelihood(candidate)
  
  lr1=logproposal(candidate,current)
  
  lr2=logproposal(current,candidate)
  
  lpcurrent=logprior(current)
  
  llcurrent=loglikelihood(current)
  
  logalpha=lpcandidate+llcandidate+lr1-lpcurrent-llcurrent-lr2
  
### 
###  Take the minimum of the log(alpha) and 0.
###


if(!is.finite(logalpha))
  {
logalpha = -Inf

  }
  
if(logalpha>0)
    {
    logalpha=0
     } 
### 
###  Generate a U(0,1) random number and take its logarithm.
### 
  logt=log(runif(1))
  
### 
###  Accept or reject the step.
### 
  if (logt < logalpha)
    {
### 
###  Accept the step.
### 
    current=candidate
    nacc=nacc+1
### 
###  Update the MAP solution if this one is better.
### 
    if ((lpcandidate+llcandidate) > lMAP)
      {
      lMAP=lpcandidate+llcandidate
      mMAP=candidate
       }
  }
  else
    {
### 
###  Reject the step.
### 
     } #
### 
###  Record the result.
### 
  mout[,k]=current
  accrate=nacc/niter
   }


return(list(mout=mout,mMAP=mMAP,accrate=accrate))

}
