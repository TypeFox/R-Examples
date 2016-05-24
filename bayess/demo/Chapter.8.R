# CHAPTER 8: R COMMANDS

# Section 8.2

data(Laichedata)
image(-t(as.matrix(Laichedata)),xlab="",ylab="",col=grey(0:99/100))

# Section 8.2.3

gisi=isingibbs(20,100,beta=0.8)
image(1:100,1:100,t(gisi),col=gray(c(1,0)),xlab="",ylab="")

S=readline(prompt="Type  <Return>   to continue : ")

hmisi=isinghm(20,100,beta=0.3)
image(1:100,1:100,t(hmisi),col=gray(c(1,0)),xlab="",ylab="")

S=readline(prompt="Type  <Return>   to continue : ")

# Section 8.2.4

hmpotts=pottshm(20,100,0.5)
image(1:100,1:100,t(hmpotts),col=gray(c(1,0.75,0.25,0)),xlab="",ylab="")

S=readline(prompt="Type  <Return>   to continue : ")

# Section 8.3.1

Z=seq(0,2,length=21)
for (i in 1:21)
Z[i]=sumising(niter=20,numb=24,beta=Z[i])
lrcst=approxfun(seq(0,2,length=21),Z)
plot(seq(0,2,length=21),Z,xlab="",ylab="")
curve(lrcst,0,2,add=T)

S=readline(prompt="Type  <Return>   to continue : ")

betatilde=0.3
beta=0.4
Zratio=integrate(lrcst,betatilde,beta)$value
Zratio

S=readline(prompt="Type  <Return>   to continue : ")

# Section 8.3.2

data(normaldata)
normaldata=normaldata[,2]
xbar=mean(normaldata)
n=length(normaldata)
s2=(n-1)*var(normaldata)
Nsim=10^6 #simulations from the prior
indem=sample(c(0,1),Nsim,rep=TRUE)
ssigma=1/rexp(Nsim)
smu=rnorm(Nsim)*sqrt(ssigma)*(indem==1)
ss2=s2/(ssigma*rchisq(Nsim,n-1))
sobs=n*(rnorm(Nsim,smu,sqrt(ssigma/n))-xbar)^2+
ss2-1-log(ss2)
epsi=quantile(sobs,.001) #bound and selection
prob=sum(indem[sobs<epsi]==0)/(0.001*Nsim)
(1-prob)/prob
 
S=readline(prompt="Type  <Return>   to continue : ")

# Section 8.3.3

ncol=4; nrow=10; Nsim=10^2; Nmc=10^2
x0=pottshm(ncol,nit=Nmc,n=nrow,beta=0.5)
suf0=(sum(x0[-1,]==x0[-nrow,])+sum(x0[,-1]==x0[,-nrow]))
outa=dista=rep(0,Nsim)
for (tt in 1:Nsim){
 beta=outa[tt]=runif(1,max=2)
 xprop=pottshm(ncol,nit=Nmc,n=nrow,beta=beta)
 dista[tt]=abs(suf0-(sum(xprop[-1,]==xprop[-nrow,])+
 sum(xprop[,-1]==xprop[,-ncol])))
 }
betas=outa[order(dista)<=.2*Nsim]

hist(betas,col="wheat2", xlab=expression(beta))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 8.4

data(Menteith)
lm3=as.matrix(Menteith)
image(1:100,1:100,lm3,col=gray(256:1/256),xlab="",ylab="")

S=readline(prompt="Type  <Return>   to continue : ")

titus=reconstruct(20,lm3)

affect=function(u)
{
order(u)[6]
}

aff=apply(titus$xcum,1,affect)
aff=t(matrix(aff,100,100))
par(mfrow=c(2,1),mar=c(2,2,2,2))
image(1:100,1:100,lm3,col=gray(256:1/256),xlab="",ylab="")
image(1:100,1:100,aff,col=gray(6:1/6),xlab="",ylab="")

S=readline(prompt="Type  <Return>   to continue : ")

par(mfrow=c(3,2),mar=c(2,2,2,2))
plot(titus$mu[,1],type="l")
plot(titus$mu[,2],type="l")
plot(titus$mu[,3],type="l")
plot(titus$mu[,4],type="l")
plot(titus$mu[,5],type="l")
plot(titus$mu[,6],type="l")

S=readline(prompt="Type  <Return>   to continue : ")

par(mfrow=c(3,2),mar=c(2,2,2,2))
hist(titus$mu[,1],nclass=100,main="")
hist(titus$mu[,2],nclass=100,main="")
hist(titus$mu[,3],nclass=100,main="")
hist(titus$mu[,4],nclass=100,main="")
hist(titus$mu[,5],nclass=100,main="")
hist(titus$mu[,6],nclass=100,main="")

S=readline(prompt="Type  <Return>   to continue : ")

par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(titus$sigma2,type="l")
hist(titus$sigma2,nclass=50,main="")
plot(titus$beta,type="l")
hist(titus$beta,nclass=50,main="")
