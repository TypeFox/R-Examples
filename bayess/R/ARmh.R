ARmh=function(x,p=1,W=10^3){
# Metropolis-Hastings algorithm for an AR(p) model

mu=mean(x)
sig2=var(x)/100

lambdareal=2*runif(p)-1
preal=p
pcomp=0
lambdacomp=0
llo=ARllog(p,x,pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,mu=mu,sig2=sig2)

psis=matrix(0,ncol=p,nrow=W)
mus=rep(0,W)
sigs=rep(0,W)
ncomp=rep(0,W)
llik=rep(0,W)

for (m in 1:W)
{

if (runif(1)<.1)
{
ind=sample(1:p,1)
if (ind<=pcomp)
{ 
ind=ind-(ind%%2==0) 
ppropreal=preal
ppropcomp=pcomp
lambpropreal=lambdareal
lambpropcomp=lambdacomp
lambpropcomp[ind]=lambdacomp[ind]+.05*rnorm(1)
lambpropcomp[ind+1]=lambdacomp[ind+1]+.05*rnorm(1)
}
else
{
ppropreal=preal
ppropcomp=pcomp
lambpropreal=lambdareal
lambpropcomp=lambdacomp
lambpropreal[ind-pcomp]=lambdareal[ind-pcomp]+.05*rnorm(1) 
}

lloprop=ARllog(p,x,pr=ppropreal,pc=ppropcomp,lr=lambpropreal,lc=lambpropcomp,mu=mu,sig2=sig2)

if (log(runif(1))<lloprop$ll-llo$ll)
{
llo=lloprop
preal=ppropreal
pcomp=ppropcomp
lambdacomp=lambpropcomp
lambdareal=lambpropreal
}
}

if (preal<2)
{
ppropreal=preal+2
ppropcomp=pcomp-2
ind=sample(1:pcomp,1)
ind=ind-(ind%%2==0) 

if (preal==0) lambpropreal=2*runif(2)-1
else lambpropreal=c(lambdareal,lambdacomp[ind]+0.01*rnorm(2))
      
lambpropcomp=lambdacomp[((1:pcomp)!=ind)&((1:pcomp)!=(ind+1))]
coef=(1/2) * (pi/4)
}

if (pcomp==0)
{
ppropreal=p-2
ppropcomp=2
ind=sample(1:p,2) 
lambpropreal=lambdareal[((1:p)!=ind[1])&((1:p)!=ind[2])] 
lambpropcomp=c(mean(lambdareal[ind])+.01*rnorm(1),.01*rnorm(1))
coef=(1/2) * (4/pi)
}

if ((preal>1)&&(pcomp>0))
{
if (runif(1)<.1)
{
ppropcomp=pcomp-2
ppropreal=preal+2
ind=sample(1:pcomp,1)
ind=ind-(ind%%2==0) 
if (ppropcomp>0) lambpropcomp=lambdacomp[((1:pcomp)!=ind)&((1:pcomp)!=(ind+1))]
else lambpropcomp=0
lambpropreal=c(lambdareal,lambdacomp[ind]+0.05*rnorm(2))
coef=9*(1+(preal<2)) * (pi/4)
}
else
{ 
ppropreal=preal-2
ppropcomp=pcomp+2
ind=sample(1:preal,2)
                                                                                                                   
if (ppropreal>0) lambpropreal=lambdareal[((1:preal)!=ind[1])&((1:preal)!=ind[2])]
else lambpropreal=0	

lambpropcomp=c(lambdacomp,mean(lambdareal[ind])+.05*rnorm(1),.05*rnorm(1))
coef=(4/pi) * (1+(ppropcomp<p-1))/9.
}
}

lloprop=ARllog(p,x,pr=ppropreal,pc=ppropcomp,lr=lambpropreal,lc=lambpropcomp,mu=mu,sig2=sig2)

if (log(runif(1))<log(coef)+lloprop$ll-llo$ll)
{
llo=lloprop
preal=ppropreal
pcomp=ppropcomp
lambdacomp=lambpropcomp
lambdareal=lambpropreal
}
psis[m,]=llo$ps[2:(p+1)]

if (runif(1)<.1)
{ 
muold=mu
mu=rnorm(1,mean=muold,sd=sqrt(2*sig2))
lloprop=ARllog(p,x,pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=llo$ps,mu=mu,sig2=sig2)
if (log(runif(1))>lloprop$ll-llo$ll) mu=muold
else llo=lloprop

sig2old=sig2
sig2=exp(rnorm(1,mean=log(sig2),sd=sqrt(2)*sig2))
lloprop=ARllog(p,x,pr=preal,pc=pcomp,lr=lambdareal,lc=lambdacomp,compsi=F,pepsi=llo$ps,mu=mu,sig2=sig2)
if (log(runif(1))>lloprop$ll-llo$ll) sig2=sig2old
else llo=lloprop
}

psis[m,]=llo$ps[2:(p+1)]
mus[m]=mu
sigs[m]=sig2
llik[m]=llo$ll
ncomp[m]=pcomp
}

list(psis=psis,mus=mus,sigs=sigs,llik=llik,ncomp=ncomp)
}
