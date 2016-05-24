MAmh=function(x,p=1,W=10^3){
# Metropolis-Hastings estimation of the parameters of an MA(p) model

mu=mean(x)
sig2=var(x)/100
varef=sig2
eps=rnorm(p,sd=sqrt(sig2))

lambdareal=2*runif(p)-1
preal=p
pcomp=0
lambdacomp=0
llo=MAllog(p,x,preal,pcomp,lambdareal,lambdacomp,mu,sig2,eps=eps)
Psi=llo$ps

# MCMC REVERSIBLE SCHEME

psis=matrix(0,ncol=p,nrow=W)
mus=rep(0,W)
sigs=rep(0,W)
ncomp=rep(0,W)
llik=rep(0,W)
epsrec=matrix(0,ncol=p,nrow=W)
preds=matrix(0,ncol=p,nrow=W)

for (m in 1:W)
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

lloprop=MAllog(p,x,ppropreal,ppropcomp,lambpropreal,lambpropcomp,mu,sig2,eps=eps)
if (log(runif(1))<lloprop$ll-llo$ll)
{
llo=lloprop
preal=ppropreal
pcomp=ppropcomp
lambdacomp=lambpropcomp
lambdareal=lambpropreal
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
lambpropcomp=c(mean(lambdareal[ind])+.01*rnorm(1),.01*rnorm(1)) #rho*c(cos(theta),sin(theta))
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

lloprop=MAllog(p,x,ppropreal,ppropcomp,lambpropreal,lambpropcomp,mu,sig2,eps=eps)
if (log(runif(1))<log(coef)+lloprop$ll-llo$ll)
{
llo=lloprop
preal=ppropreal
pcomp=ppropcomp
lambdacomp=lambpropcomp
lambdareal=lambpropreal
}
psis[m,]=llo$ps
muold=mu
mu=rnorm(1,mean=muold,sd=sqrt(5*sig2))
lloprop=MAllog(p,x,preal,pcomp,lambdareal,lambdacomp,mu,sig2,compsi=F,pepsi=llo$ps,eps=eps)
if (log(runif(1))>lloprop$ll-llo$ll) mu=muold
else llo=lloprop

if (runif(1)<.3)
{
Psi=llo$ps
heps=rep(0,2*p)
heps[1:p]=eps  
for (j in (p+1):(2*p)) heps[j]=x[j]-sum(Psi*heps[(j-p):(j-1)])
varheps=(2*p)*var(heps)
sig2old=sig2
sig2=varheps/rgamma(1,2*p)
difdens=dgamma(sig2/varheps,2*p,log=T)-dgamma(sig2old/varheps,2*p,log=T)+log(sig2)-log(sig2old)
lloprop=MAllog(p,x,preal,pcomp,lambdareal,lambdacomp,mu,sig2,compsi=F,pepsi=llo$ps,eps=eps)
if (log(runif(1))>lloprop$ll-llo$ll-difdens) sig2=sig2old
else llo=lloprop
}

if (runif(1)<.3)
{
sig2old=sig2
sig2=exp(rnorm(1,mean=log(sig2),sd=sqrt(.1*varef)))
lloprop=MAllog(p,x,preal,pcomp,lambdareal,lambdacomp,mu,sig2,compsi=F,pepsi=llo$ps,eps=eps)
difdens=log(sig2)-log(sig2old)
if (log(runif(1))>lloprop$ll-llo$ll-difdens) sig2=sig2old
else llo=lloprop
}
if (runif(1)<.3)
{
sig2old=sig2
thismean=varef/(1+sum(Psi^2))
sig2=exp(rnorm(1,mean=thismean,sd=sqrt(varef)))
difdens=dnorm(log(sig2),mean=thismean,sd=sqrt(.5*varef),log=T)-dnorm(log(sig2old),mean=thismean,sd=sqrt(.5*varef),log=T)
lloprop=MAllog(p,x,preal,pcomp,lambdareal,lambdacomp,mu,sig2,compsi=F,pepsi=llo$ps,eps=eps)
if (log(runif(1))>lloprop$ll-llo$ll-difdens) sig2=sig2old
else llo=lloprop
}

Psi=llo$ps

if (runif(1)<.5)
{
heps=rep(0,2*p+1)
keps=rep(0,p)
for (i in 1:p)
{
x=x-mu
heps[1:p]=eps
for (j in (p+1):(2*p+1)) heps[j]=x[j]+sum(rev(Psi)*heps[(j-p):(j-1)])
heps[i]=0
for (j in 1:(p-i+1)) keps[j]=x[j]+sum(rev(Psi)*heps[j:(j+p-1)])
x=x+mu
epsvar=1/sum(c(1,Psi[i:p]^2))
epsmean=sum(Psi[i:p]*keps[1:(p-i+1)])*epsvar
epsmean=epsmean/epsvar
epsvar=sig2*epsvar
propeps=rnorm(1,mean=epsmean,sd=sqrt(epsvar))
epspr=eps
epspr[i]=propeps
lloprop=MAllog(p,x,preal,pcomp,lambdareal,lambdacomp,mu,sig2,compsi=F,pepsi=Psi,eps=epspr)
propsal1=dnorm(propeps,mean=epsmean,sd=sqrt(epsvar),log=T)
x=x-mu
heps[i]=propeps
for (j in (p+1):(2*p+1)) heps[j]=x[j]+sum(rev(Psi)*heps[(j-p):(j-1)])
heps[i]=0
for (j in 1:(p-i+1)) keps[j]=x[j]+sum(rev(Psi)*heps[j:(j+p-1)])
x=x+mu
epsvar=1/sum(c(1,Psi[i:p]^2))
epsmean=sum(Psi[i:p]*keps[1:(p-i+1)])
epsmean=epsmean*epsvar
epsvar=sig2*epsvar 
propsal0=dnorm(eps[i],mean=epsmean,sd=sqrt(epsvar),log=T)
if (log(runif(1))<lloprop$ll-llo$ll-propsal1+propsal0)
{
eps[i]=propeps;
llo=lloprop
}
}
}
if (runif(1)<.5)
{
heps=rep(0,2*p+1)
keps=rep(0,p)
for (i in 1:p)
{ 
propeps = rnorm(1,mean=eps[i],sd=0.1*sqrt(sig2))
heps[1:p]=eps
heps[i]=propeps
keps=heps[1:p]
allop=MAllog(p,x,preal,pcomp,lambdareal,lambdacomp,mu,sig2,compsi=F,pepsi=Psi,eps=keps)
if (log(runif(1))<allop$ll-llo$ll)
{ 
eps[i]=propeps;
llo=allop
}
}
}
psis[m,]=llo$ps
mus[m]=mu
sigs[m]=sig2
llik[m]=llo$ll
ncomp[m]=pcomp
epsrec[m,]=eps
heps=rep(0,length(x)+2*p)
heps[1:p]=eps
hatx=heps
hatx[1:length(x)]=x-mu
for (j in (p+1):(p+length(x))) heps[j]=x[j-p]+sum(rev(Psi)*heps[(j-p):(j-1)])
for (j in (p+length(x)+1):(p+length(x)+p)) hatx[j-p]=-sum(rev(Psi)*heps[(j-p):(j-1)])
preds[m,]=hatx[(length(x)+1):(length(x)+p)]+mu
}

list(psis=psis,mus=mus,sigs=sigs,llik=llik,ncomp=ncomp,epsrec=epsrec)
}
