MAllog=function(p,dat,pr,pc,lr,lc,mu,sig2,compsi=T,pepsi=rep(0,p),eps=rnorm(p))
{

# LIKELIHOOD REPRESENTATION

Psi=matrix(0,ncol=p,nrow=p+1)

if (compsi)
{

Psi[1,]=1

if (pr>0)
{
Psi[2,1]=-lr[1]
if (pr>1)
{
for (i in 2:pr) Psi[2:(i+1),i]=Psi[2:(i+1),i-1]-lr[i]*Psi[1:i,i-1]
}
}
if (pc>0)
{
if (pr>0)
{
Psi[2,pr+2]=-2*lc[1]+Psi[2,pr]
Psi[3:(pr+3),pr+2]=(lc[1]^2+lc[2]^2)*Psi[1:(pr+1),pr]-2*lc[1]*Psi[2:(pr+2),pr]+Psi[3:(pr+3),pr]
}
else
{
Psi[2,2]=-2*lc[1];
Psi[3,2]=(lc[1]^2+lc[2]^2);
}
if (pc>2)
{
for (i in seq(4,pc,2))
{
pri=pr+i
prim=pri-2
Psi[2,pri]=-2*lc[i-1]+Psi[2,prim]
Psi[3:(pri+1),pri]=(lc[i-1]^2+lc[i]^2)*Psi[1:(pri-1),prim]-2*lc[i-1]*Psi[2:pri,prim]+Psi[3:(pri+1),prim]
}
}
}
Psi=Psi[2:(p+1),p]
}
else
{
Psi=pepsi
}

# LIKELIHOOD CALCULATION

dat=dat-mu
heps=rep(0,length(dat)+p)
heps[1:p]=eps
for (i in 1:length(dat))
heps[p+i]=dat[i]+sum(rev(Psi)*heps[i:(p+i-1)])
loglike=-((sum(heps^2)/sig2)+(length(dat)+p)*log(sig2))/2

if (loglike==Inf || is.nan(loglike) || loglike==-Inf) loglike=-1e300

list(ll=loglike,ps=Psi)
}
