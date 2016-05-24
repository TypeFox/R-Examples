# Chapter 6 R commands

library(combinat)

data(datha)
datha=as.matrix(datha)
hist(datha,nclas=200,xlab="",xlim=c(min(datha),max(datha)),ylab="",prob=TRUE,main="")

S=readline(prompt="Type  <Return>   to continue : ")

# Section 6.3

n=length(datha)
meand=mean(datha)
vard=var(datha)
sdd=sqrt(vard)

dat=plotmix()$sample

S=readline(prompt="Type  <Return>   to continue : ")

omega=function(z,x,p)
{
n=length(x)
n1=sum(z==1)
n2=n-n1
if (n1==0) xbar1=0 else xbar1=sum((z==1)*x)/n1
if (n2==0) xbar2=0 else xbar2=sum((z==2)*x)/n2
ss1=sum((z==1)*(x-xbar1)^2)
ss2=sum((z==2)*(x-xbar2)^2)
sqrt((n1+.25)*(n2+.25))*p^n1*(1-p)^n2*
exp(-((n1+.25)*ss1+(n2+.25)*ss2)/2)*
exp(-(n1*xbar1^2+n2*xbar2)/8)
}
omega(z=sample(1:2,4,rep=TRUE),x=plotmix(n=4,plot=FALSE)$samp,p=.8)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 6.4

dat=plotmix()$sample
simu=gibbsmean(0.7,dat)
points(simu,pch=".")

S=readline(prompt="Type  <Return>   to continue : ")

mix=list(k=3,mu=mean(datha),sig=var(datha))
simu=gibbsnorm(100,datha,mix)
hist(datha,prob=T,main="",xlab="",ylab="",nclass=100,col="wheat")
x=y=seq(min(datha),max(datha),length=150)
yy=matrix(0,ncol=150,nrow=100)
for (i in 1:150)
{
yy[,i]=apply(simu$p*dnorm(x[i],mean=simu$mu,sd=sqrt(simu$sig)),1,sum)
y[i]=mean(yy[,i])
}
for (t in 51:100) lines(x,yy[t,],col="gold")
lines(x,y,lwd=2.3,col="sienna2")

S=readline(prompt="Type  <Return>   to continue : ")

dat=plotmix()$sample
simean=hmmeantemp(dat,10000,var=1)
points(simean,pch=19,cex=.4)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 6.5

# recycling former simu
par(mfrow=c(3,2),mar=c(4,4,1,1.5))
plot(simu$mu[,1],ylim=range(simu$mu),ylab=expression(mu[i]),xlab="n",type="l",col="gold4")
lines(simu$mu[,2],col="sienna3")
lines(simu$mu[,3],col="steelblue")
plot(simu$mu[,2],simu$p[,2],col="sienna3",
xlim=range(simu$mu),ylim=range(simu$p), 
xlab=expression(mu[i]),ylab=expression(p[i]))
points(simu$mu[,3],simu$p[,3],col="steelblue")
points(simu$mu[,1],simu$p[,1],col="gold4")
plot(simu$p[,1],ylim=range(simu$p),ylab=expression(p[i]),xlab="n",type="l",col="gold4")
lines(simu$p[,2],col="sienna3")
lines(simu$p[,3],col="steelblue")
plot(simu$p[,2],simu$sig[,2],col="sienna3",
xlim=range(simu$p),ylim=range(simu$sig), 
xlab=expression(p[i]),ylab=expression(sigma[i]^2))
points(simu$p[,3],simu$sig[,3],col="steelblue")
points(simu$p[,1],simu$sig[,1],col="gold4")
plot(simu$sig[,1],ylim=range(simu$sig),ylab=expression(sigma[i]^2),xlab="n",type="l",col="gold4")
lines(simu$sig[,2],col="sienna3")
lines(simu$sig[,3],col="steelblue")
plot(simu$sig[,2],simu$mu[,2],col="sienna3",
xlim=range(simu$sig),ylim=range(simu$mu), 
xlab=expression(sigma[i]^2),ylab=expression(mu[i]))
points(simu$sig[,3],simu$mu[,3],col="steelblue")
points(simu$sig[,1],simu$mu[,1],col="gold4")

S=readline(prompt="Type  <Return>   to continue : ")

indimap=order(simu$lopost,decreasing=TRUE)[1]
map=list(mu=simu$mu[indimap,],
sig=simu$sig[indimap,],
p=simu$p[indimap,])
lili=alloc=matrix(0,length(datha),3)
for (t in 1:length(datha))
{
lili[t,]=map$p*dnorm(datha[t],mean=map$mu,
sd=sqrt(map$sig))
lili[t,]=lili[t,]/sum(lili[t,])
}

ormu=orsig=orp=matrix(0,ncol=3,nrow=1000)
perma=permn(3)
for (t in 1:1000)
{
entropies=rep(0,factorial(3))
for (j in 1:n)
{
alloc[j,]=simu$p[t,]*dnorm(datha[j],mean=simu$mu[t,],
sd=sqrt(simu$sig[t,]))
alloc[j,]=alloc[j,]/sum(alloc[j,])
for (i in 1:factorial(3))
entropies[i]=entropies[i]+
sum(lili[j,]*log(alloc[j,perma[[i]]]))
}
best=order(entropies,decreasing=TRUE)[1]
ormu[t,]=simu$mu[t,perma[[best]]]
orsig[t,]=simu$sig[t,perma[[best]]]
orp[t,]=simu$p[t,perma[[best]]]
}


# Section 6.7

mu1=2.5;mu2=0;p=.7;n=500;nl=50
pbar=1-p
u=runif(n)
sampl=rnorm(n)+(u<=p)*mu1+(u>p)*mu2 
mu1=mu2=seq(min(sampl),max(sampl),.1) 
mo1=mu1%*%t(rep(1,length(mu2))) 
mo2=rep(1,length(mu2))%*%t(mu2) 
ca1=-0.5*mo1*mo1 
ca2=-0.5*mo2*mo2 
like=0*mo1 
for (i in 1:n) 
like=like+log(p*exp(ca1+sampl[i]*mo1)+ 
pbar*exp(ca2+sampl[i]*mo2)) 
like=like+.1*(ca1+ca2) 

par(mfrow=c(1,3),mar=c(4,4,1,1))
for (i in 1:3){
 image(mu1,mu2,like,xlab=expression(mu[1]), 
 ylab=expression(mu[2]),col=heat.colors(250)) 
 contour(mu1,mu2,like,levels=seq(min(like),max(like),length=nl),
 add=TRUE,drawlabels=FALSE) 
 simu=hmmeantemp(dat,10000,var=0.1,alpha=1/(10)^(i-1))
 points(simu,pch=19,cex=.4)
 }

S=readline(prompt="Type  <Return>   to continue : ")

dat=plotmix()$sample
simu=hmmeantemp(dat,10000)
points(simu,pch=19,cex=.4,col="sienna")
simu=hmmeantemp(dat,10000,alpha=0.1)
points(simu,pch=19,cex=.4,col="steelblue4")
simu=hmmeantemp(dat,10000,alpha=0.01)
points(simu,pch=19,cex=.4)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 6.8

########## k=2

k=2
mu=rnorm(k,mean=meand,sd=sdd/k)
sig=1/rgamma(k,shape=10,scale=vard)
p=rdirichlet(par=rep(1,k))
mix=list(k=k,p=p,mu=mu,sig=sig)
resim2=gibbs(200,datha,mix)
lulu=order(resim2$lolik)[200]
lnum1=resim2$lolik[lulu]
lnum2=sum(dnorm(resim2$mu[lulu,],mean=meand,sd=resim2$sig[lulu,],log=TRUE)+
dgamma(1/resim2$sig[lulu,],10,var(datha),log=TRUE)-2*log(resim2$sig[lulu,]))+
sum((rep(0.5,2)-1)*log(resim2$p[lulu,]))+
lgamma(sum(rep(0.5,2)))-sum(lgamma(rep(0.5,2)))
lchibapprox2=lnum1+lnum2-log(resim2$deno)
lchibapprox2

S=readline(prompt="Type  <Return>   to continue : ")

########## k=3

k=3
mu=rnorm(k,mean=meand,sd=sdd/k)
sig=1/rgamma(k,shape=10,scale=vard)
p=rdirichlet(par=rep(1,k))
mix=list(k=k,p=p,mu=mu,sig=sig)
resim3=gibbs(200,datha,mix)
lulu=order(resim3$lolik)[200]
lnum1=resim3$lolik[lulu]
lnum2=sum(dnorm(resim3$mu[lulu,],mean=meand,sd=resim3$sig[lulu,],log=TRUE)+
dgamma(resim3$sig[lulu,],10,vard,log=TRUE)-2*log(resim3$sig[lulu,]))+
sum((rep(0.5,mix$k)-1)*log(resim3$p[lulu,]))+
lgamma(sum(rep(0.5,mix$k)))-sum(lgamma(rep(0.5,mix$k)))
lchibapprox3=lnum1+lnum2-log(resim3$deno)

S=readline(prompt="Type  <Return>   to continue : ")

########## k=4

k=4
mu=rnorm(k,mean=meand,sd=sdd/k)
sig=1/rgamma(k,shape=10,scale=vard)
p=rdirichlet(par=rep(1,k))
mix=list(k=k,p=p,mu=mu,sig=sig)
resim4=gibbs(200,datha,mix)
lulu=order(resim4$lolik)[200]
lnum1=resim4$lolik[lulu]
lnum2=sum(dnorm(resim4$mu[lulu,],mean=meand,sd=resim4$sig[lulu,],log=TRUE)+
dgamma(resim4$sig[lulu,],10,vard,log=TRUE)-2*log(resim4$sig[lulu,]))+
sum((rep(0.5,mix$k)-1)*log(resim4$p[lulu,]))+
lgamma(sum(rep(0.5,mix$k)))-sum(lgamma(rep(0.5,mix$k)))
lchibapprox4=lnum1+lnum2-log(resim4$deno)

S=readline(prompt="Type  <Return>   to continue : ")

########## k=5

k=5
mu=rnorm(k,mean=meand,sd=sdd/k)
sig=1/rgamma(k,shape=10,scale=vard)
p=rdirichlet(par=rep(1,k))
mix=list(k=k,p=p,mu=mu,sig=sig)
resim5=gibbs(200,datha,mix)
lulu=order(resim5$lolik)[200]
lnum1=resim5$lolik[lulu]
lnum2=sum(dnorm(resim5$mu[lulu,],mean=meand,sd=resim5$sig[lulu,],log=TRUE)+
dgamma(resim5$sig[lulu,],10,vard,log=TRUE)-2*log(resim5$sig[lulu,]))+
sum((rep(0.5,mix$k)-1)*log(resim5$p[lulu,]))+
lgamma(sum(rep(0.5,mix$k)))-sum(lgamma(rep(0.5,mix$k)))
lchibapprox5=lnum1+lnum2-log(resim5$deno)

S=readline(prompt="Type  <Return>   to continue : ")

########## k=6

k=6
mu=rnorm(k,mean=meand,sd=sdd/k)
sig=1/rgamma(k,shape=10,scale=vard)
p=rdirichlet(par=rep(1,k))
mix=list(k=k,p=p,mu=mu,sig=sig)
resim6=gibbs(200,datha,mix)
lulu=order(resim6$lolik)[200]
lnum1=resim6$lolik[lulu]
lnum2=sum(dnorm(resim6$mu[lulu,],mean=meand,sd=resim6$sig[lulu,],log=TRUE)+
dgamma(resim6$sig[lulu,],10,vard,log=TRUE)-2*log(resim6$sig[lulu,]))+
sum((rep(0.5,mix$k)-1)*log(resim6$p[lulu,]))+
lgamma(sum(rep(0.5,mix$k)))-sum(lgamma(rep(0.5,mix$k)))
lchibapprox6=lnum1+lnum2-log(resim6$deno)
