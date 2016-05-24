TAR.simu<-function(nob,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,thres,lagp1,lagp2){

## Simulate nmax observations and discard first
## 1000 observations.
nmax<-nob+1000                ## No. of simulated observations
p<-max(p1,p2)+1                         ## Set initial period
y<-rep(0,nmax)                ## Generate nmax zeros to y
at<-rnorm(nmax,mean=0,sd=1)   ## Simulate nmax normal random variables

for (i in p:nmax){               ## Generate data y recursively from time p to nmax
if (y[i-lagd] <= thres)          ## Determine the location of the lagged y
                                 ## Simulate data from Regime 1
{y[i]=sum(ph.1[1],y[i-lagp1]*ph.1[2:(p1+1)],at[i]*sqrt(sig.1))}
else
                                 ## Simulate data from Regime 2
{y[i]=sum(ph.2[1],y[i-lagp2]*ph.2[2:(p2+1)],at[i]*sqrt(sig.2))}
}

yt<-y[1001:nmax]             ## Discard first 1000 observations and
                             ## save the observations from 1001 to nmax.
return(yt)
}

