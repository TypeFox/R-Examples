# Chapter 2 R commands

# Section 2.1

data(normaldata)
shift=normaldata[,2]
hist(shift,nclass=10,col="steelblue",prob=TRUE,main="")

S=readline(prompt="Type  <Return>   to continue : ")

qqnorm((shift-mean(shift))/sd(shift),pch=19,col="gold2")
abline(a=0,b=1,lty=2,col="indianred",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.2.1

n=length(shift)
mmu=sum(shift)/(n+1);mmu
vmu=0.75^2/(n+1);vmu

S=readline(prompt="Type  <Return>   to continue : ")

mtmu=sum(shift)/(n+1);mtmu
stmu=(1+(n-1)*var(shift))/((n+2)*(n+1));stmu

S=readline(prompt="Type  <Return>   to continue : ")

library(mnormt)
curve(dmt(x,mean=mmu,S=stmu,df=n+2),col="chocolate2",lwd=2,
xlab="x",ylab="",xlim=c(-.5,.5))
curve(dnorm(x,mean=mmu,sd=sqrt(vmu)),col="steelblue2",
lwd=2,add=TRUE,lty=2)

S=readline(prompt="Type  <Return>   to continue : ")

digmma=function(x,shape,scale){
   dgamma(1/x,shape,scale)/x^2}
curve(digmma(x,shape=33,scale=(1+(n+1)*var(shift))/2),xlim=c(0,.2),lwd=2)   
pgamma(1/(.75)^2,shape=33,scale=(1+(n+1)*var(shift))/2)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.2.5

qt(.975,df=n)*sqrt((n-1)*var(shift)/n^2)
-qt(.975,df=n)*sqrt((n-1)*var(shift)/n^2)+mean(shift)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.3.2

BaFa=function(z,rat){ 
   sqrt(1/(1+rat))*exp(z^2/(2*(1+1/rat)))}
BaFa(mean(shift),1)
BaFa(mean(shift),10)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.3.3

ratio=n*mean(shift)^2/((n-1)*var(shift))
((1+ratio)/(1+ratio/(n+1)))^(n/2)/sqrt(n+1)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.4.1

illing=normaldata
xsam=illing[illing[,1]==5,2]
xbar=mean(xsam);xbar
ysam=illing[illing[,1]==6,2]
ybar=mean(ysam);ysam
Ssquar=9*(var(xsam)+var(ysam))/10;Ssquar

S=readline(prompt="Type  <Return>   to continue : ")

Nsim=10^4
tau=0.75
xis=rnorm(Nsim,sd=tau)
BaFa=mean(((2*xis+xbar-ybar)^2+2*Ssquar)^(-8.5))/((xbar-ybar)^2+2*Ssquar)^(-8.5)
xis=matrix(rnorm(500*10^3,sd=tau),nrow=500);BaFa

S=readline(prompt="Type  <Return>   to continue : ")

BF=((2*xis+xbar-ybar)^2+2*Ssquar)^(-8.5)/((xbar-ybar)^2+2*Ssquar)^(-8.5)
estims=apply(BF,1,mean) 
hist(estims,nclass=84,prob=T,col="wheat2",main="",xlab="Bayes Factor estimates") 
curve(dnorm(x,mean=mean(estims),sd=sd(estims)),col="steelblue2",add=TRUE)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.4.2

Nobs=10
obs=rcauchy(Nobs)
Nsim=250
Nmc=500
sampl=matrix(rnorm(Nsim*Nmc),nrow=1000) # normal samples
raga=riga=matrix(0,nrow=50,ncol=2) # ranges
mu=0
for (j in 1:50){
 prod=1/dnorm(sampl-mu) # importance sampling
 for (i in 1:Nobs)
 prod=prod*dt(obs[i]-sampl,1)
 esti=apply(sampl*prod,2,sum)/apply(prod,2,sum)
 raga[j,]=range(esti)
 riga[j,]=c(quantile(esti,.025),quantile(esti,.975))
 sampl=sampl+0.1
 mu=mu+0.1
 }
mus=seq(0,4.9,by=0.1)
plot(mus,0*mus,col="white",xlab=expression(mu),
 ylab=expression(hat(theta)),ylim=range(raga))
polygon(c(mus,rev(mus)),c(raga[,1],rev(raga[,2])),col="grey50")
polygon(c(mus,rev(mus)),c(riga[,1],rev(riga[,2])),col="pink3")

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.4.3

((1+ratio)/(1+ratio/(n+1)))^(-n/2)*sqrt(n+1)

S=readline(prompt="Type  <Return>   to continue : ")

n=64
xbar=mean(shift)
sqar=(n-1)*var(shift)
Nmc=10^7
# Simulation from model M2:
sigma2=1/rgamma(Nmc,shape=n/2,rate=(n*xbar^2/(n+1)+sqar)/2)
mu2=rnorm(Nmc,n*xbar/(n+1),sd=sqrt(sigma2/(n+1)))
# Simulation from model M1:
sigma1=1/rgamma(Nmc,shape=n/2,rate=(n*xbar^2+sqar)/2)
muhat=mean(mu2)
tauat=sd(mu2)
mu1=rnorm(Nmc,mean=muhat,sd=tauat)
#tilde functions
tildepi1=function(sigma,mu){
 exp(-.5*((n*xbar^2+sqar)/sigma+(n+2)*log(sigma))+
 dnorm(mu,muhat,tauat,log=T))
 }
 tildepi2=function(sigma,mu){
 exp(-.5*((n*(xbar-mu)^2+sqar+mu^2)/sigma+(n+3)*log(sigma)+
 log(2*pi)))}
#Bayes Factor loop
K=diff=1
rationum=tildepi2(sigma1,mu1)/tildepi1(sigma1,mu1)
ratioden=tildepi1(sigma2,mu2)/tildepi2(sigma2,mu2)
while (diff>0.01*K){
 BF=mean(1/(1+K*rationum))/mean(1/(K+ratioden))
 diff=abs(K-BF)
 K=BF}
BF

S=readline(prompt="Type  <Return>   to continue : ")

sigma1=1/rgamma(Nmc,shape=n/2,rate=(n*xbar^2+sqar)/2)
sihat=mean(log(sigma1))
tahat=sd(log(sigma1))
sigma1b=exp(rnorm(Nmc,sihat,tahat))
#tilde function
tildepi1=function(sigma){
  exp(-.5*((n*xbar^2+sqar)/sigma+(n+2)*log(sigma)))
  }
K=diff=1
rnum=dnorm(log(sigma1b),sihat,tahat)/(sigma1b*tildepi1(sigma1b))
rden=sigma1*tildepi1(sigma1)/dnorm(log(sigma1),sihat,tahat)
while (diff>0.01*K){
BF=mean(1/(1+K*rnum))/mean(1/(K+rden))
diff=abs(K-BF)
K=BF}
m1=BF

sigma2=1/rgamma(Nmc,shape=n/2,rate=(n*xbar^2/(n+1)+sqar)/2)
mu2=rnorm(Nmc,n*xbar/(n+1),sd=sqrt(sigma2/(n+1)))
temean=c(mean(mu2),mean(log(sigma2)))
tevar=cov.wt(cbind(mu2,log(sigma2)))$cov
te2b=rmnorm(Nmc,mean=temean,tevar)
mu2b=te2b[,1]
sigma2b=exp(te2b[,2])
tildepi2=function(sigma){
   exp(-.5*((n*xbar^2+sqar)/sigma+(n+2)*log(sigma)))
   }
K=diff=1
rnum=dnorm(log(sigma2b),sihat,tahat)/(sigma2b*tildepi2(sigma2b))
rden=sigma2*tildepi2(sigma2)/dnorm(log(sigma2),sihat,tahat)
while (diff>0.01*K){
BF=mean(1/(1+K*rnum))/mean(1/(K+rden))
diff=abs(K-BF)
K=BF}
m2=BF

m1/m2

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.5

n=length(shift)
outl=rep(0,n)
for (i in 1:n){
 lomean=-mean(shift[-i])
 losd=sd(shift[-i])*sqrt((n-2)/n)
 outl[i]=pt((shift[i]-lomean)/losd,df=n-1)
 }

plot(c(0,1),c(0,1),lwd=2,ylab="Predictive",xlab="Uniform",type="l")
points((1:n)/(n+1),sort(outl),pch=19,col="steelblue3")
points((1:n)/(n+1),sort(runif(n)),pch=19,col="tomato")
