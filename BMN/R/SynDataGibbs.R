############################################
###
### Gibbs sampler for synthetic data
###
############################################

GibbsIteration = function(x, Theta)
{
    n=length(x)
    for(i in 1:n)
    {
        x[i]=1
        foo=exp(sum(Theta[i,]*x))
        p=foo/(foo+1)
        x[i]=rbinom(1,1,prob=p)
    }
    return(x)
}


nextSample = function(x,Theta, iters)
{
    for(i in 1:iters)
    {
        x=GibbsIteration(x,Theta)
    }
    return(x)
}

BMNSamples=function(Theta, numSamples, burnIn, skip)
{
    p=dim(Theta)[1]
    x=matrix(numeric(p*numSamples), nrow=numSamples)
    x[1,] = nextSample(rep(0,p),Theta,burnIn)
    if(numSamples>1)
    {
        for(i in 2:numSamples)
        {
            x[i,] = nextSample(x[i-1,],Theta,skip)
        }
    }
    return(x)
}

###########################################################
###
### write a function that finds the second moment matrix using 
### gibbs sampling
###
###########################################################

gibbsSecMoment=function(Theta, numSamples=10^5)
{
    X=BMNSamples(Theta, numSamples, 1, 1)
    return(t(X) %*% X / numSamples)
}

