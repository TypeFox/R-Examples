challenge=function(Nsim=10^4){
#Slice sampler for the Challenger data

rtrun=function(mu=0,sigma=1,minn=-Inf,maxx=Inf){
#Simulation from a truncated normal
    qnorm(pnorm(minn,mean=mu,sd=sigma)+(pnorm(maxx,mean=mu,sd=sigma)-pnorm(minn,mean=mu,sd=sigma))*runif(1))
    }

#data(challenger)
x=challenger[,2]
y=challenger[,1]
n=length(x)

a=b=rep(0,Nsim)
#starts with mle:
a[1]=glm(challenger[,1]~challenger[,2])$coef[1]
b[1]=glm(challenger[,1]~challenger[,2])$coef[2]

#normal priors
meana=0;sigmaa=5
meanb=0;sigmab=5/sd(x)

for (t in 2:Nsim){

  #auxiliary uniforms
  uni=runif(n)*exp(y*(a[t-1]+b[t-1]*x))/(1+exp(a[t-1]+b[t-1]*x))

  mina=max(log(uni[y==1]/(1-uni[y==1]))-b[t-1]*x[y==1])
  maxa=min(-log(uni[y==0]/(1-uni[y==0]))-b[t-1]*x[y==0])
  a[t]=rtrun(0,sigmaa,mina,maxa)

  #auxiliary uniforms
  uni=runif(n)*exp(y*(a[t]+b[t-1]*x))/(1+exp(a[t]+b[t-1]*x))

  minb=max((log(uni[y==1]/(1-uni[y==1]))-a[t])/x[y==1])
  maxb=min((-log(uni[y==0]/(1-uni[y==0]))-a[t])/x[y==0])
  b[t]=rtrun(0,sigmab,minb,maxb)
  }
  list(a=a,b=b)
  }
