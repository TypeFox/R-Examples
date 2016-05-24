# R commands for Chapter 4

library(mvtnorm)

# Section 4.1.1

data(bank)
bank=as.matrix(bank)
y=bank[,5]
X=bank[,1:4]

par(mfrow=c(1,2))
plot(y~X[,4],xlab="Bottom margin width (mm)",ylab="Status")
boxplot(X[,4]~y,ylab="Bottom margin width (mm)",xlab="Status")

S=readline(prompt="Type  <Return>   to continue : ")

# Section 4.2.4

hm1=function(n,x0,sigma2)
#### Gaussian random walk MH sampler for a standard normal distribution
{
x=rep(0,n)
x[1]=x0
for (i in 2:n)
{
y=rnorm(1,x[i-1],sqrt(sigma2))
if (runif(1)<=exp(-0.5*(y^2-x[i-1]^2))) x[i]=y else x[i]=x[i-1]
}
x
}

test=hm1(10000,1,0.001)
par(mfrow=c(3,1))
plot(test,type="l",xlab="Iterations",ylab="MH chain")
hist(test[8001:10000],nclass=50,prob=TRUE,main="",xlab="x")
curve(dnorm,-3,3,add=TRUE)
acf(test,lag=1000,main="",ylab="Autocorrelation",ci=F)

S=readline(prompt="Type  <Return>   to continue : ")

test=hm1(10000,1,1000)
par(mfrow=c(3,1))
plot(test,type="l",xlab="Iterations",ylab="MH chain")
hist(test[8001:10000],nclass=50,prob=TRUE,main="",xlab="x")
curve(dnorm,-3,3,add=TRUE)
acf(test,lag=1000,main="",ylab="Autocorrelation",ci=F)

S=readline(prompt="Type  <Return>   to continue : ")

test=hm1(10000,1,1)
par(mfrow=c(3,1))
plot(test,type="l",xlab="Iterations",ylab="MH chain")
hist(test[8001:10000],nclass=50,prob=TRUE,main="",xlab="x")
curve(dnorm,-3,3,add=TRUE)
acf(test,lag=1000,main="",ylab="Autocorrelation",ci=F)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 4.3

# bank DATASET: GAUSSIAN RANDOM WALK MH SAMPLER UNDER FLAT PRIOR

flatprobit=hmflatprobit(10000,y,X,1)
par(mfrow=c(4,3),mar=1+c(1.5,1.5,1.5,1.5))
plot(flatprobit[,1],type="l",xlab="Iterations",ylab=expression(beta[1]))
hist(flatprobit[1001:10000,1],nclass=50,prob=TRUE,main="",xlab=expression(beta[1]))
acf(flatprobit[1001:10000,1],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(flatprobit[,2],type="l",xlab="Iterations",ylab=expression(beta[2]))
hist(flatprobit[1001:10000,2],nclass=50,prob=TRUE,main="",xlab=expression(beta[2]))
acf(flatprobit[1001:10000,2],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(flatprobit[,3],type="l",xlab="Iterations",ylab=expression(beta[3]))
hist(flatprobit[1001:10000,3],nclass=50,prob=TRUE,main="",xlab=expression(beta[3]))
acf(flatprobit[1001:10000,3],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(flatprobit[,4],type="l",xlab="Iterations",ylab=expression(beta[4]))
hist(flatprobit[1001:10000,4],nclass=50,prob=TRUE,main="",xlab=expression(beta[4]))
acf(flatprobit[1001:10000,4],lag=1000,main="",ylab="Autocorrelation",ci=F)

S=readline(prompt="Type  <Return>   to continue : ")

mean(flatprobit[1001:10000,1])
mean(flatprobit[1001:10000,2])
mean(flatprobit[1001:10000,3])
mean(flatprobit[1001:10000,4])
pnorm(-1.2193*214.9+0.9540*130.1+0.9795*129.9+1.1481*9.5)

S=readline(prompt="Type  <Return>   to continue : ")

noinfprobit=hmnoinfprobit(10000,y,X,1)
par(mfrow=c(4,3),mar=1+c(1.5,1.5,1.5,1.5))
plot(noinfprobit[,1],type="l",xlab="Iterations",ylab=expression(beta[1]))
hist(noinfprobit[1001:10000,1],nclass=50,prob=TRUE,main="",xlab=expression(beta[1]))
acf(noinfprobit[1001:10000,1],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(noinfprobit[,2],type="l",xlab="Iterations",ylab=expression(beta[2]))
hist(noinfprobit[1001:10000,2],nclass=50,prob=TRUE,main="",xlab=expression(beta[2]))
acf(noinfprobit[1001:10000,2],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(noinfprobit[,3],type="l",xlab="Iterations",ylab=expression(beta[3]))
hist(noinfprobit[1001:10000,3],nclass=50,prob=TRUE,main="",xlab=expression(beta[3]))
acf(noinfprobit[1001:10000,3],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(noinfprobit[,4],type="l",xlab="Iterations",ylab=expression(beta[4]))
hist(noinfprobit[1001:10000,4],nclass=50,prob=TRUE,main="",xlab=expression(beta[4]))
acf(noinfprobit[1001:10000,4],lag=1000,main="",ylab="Autocorrelation",ci=F)

S=readline(prompt="Type  <Return>   to continue : ")

mean(noinfprobit[1001:10000,1])
mean(noinfprobit[1001:10000,2])
mean(noinfprobit[1001:10000,3])
mean(noinfprobit[1001:10000,4])

S=readline(prompt="Type  <Return>   to continue : ")

# full model
mkprob=apply(noinfprobit,2,mean)
vkprob=var(noinfprobit)
simk=rmvnorm(10^5,mkprob,2*vkprob)
usk=probitnoinflpost(simk,y,X)-
    dmvnorm(simk,mkprob,2*vkprob,log=T)
# null model
noinfprobit0=hmnoinfprobit(10^4,y,X[,3:4],1)
mk0=apply(noinfprobit0,2,mean)
vk0=var(noinfprobit0)
simk0=rmvnorm(10^4,mk0,2*vk0)
usk0=probitnoinflpost(simk0,y,X[,3:4])-
     dmvnorm(simk0,mk0,2*vk0,log=T)
# Bayes factor
bf0probit=mean(exp(usk))/mean(exp(usk0))
bf0probit

S=readline(prompt="Type  <Return>   to continue : ")

mkprob=apply(noinfprobit,2,mean)
vkprob=var(noinfprobit)
simk=rmvnorm(10^5,mkprob,2*vkprob)
usk=probitnoinflpost(simk,y,X)-dmvnorm(simk,mkprob,2*vkprob,log=T)

noinfprobit0=hmnoinfprobit(10000,y,X[,3:4],1)
mk0=apply(noinfprobit0,2,mean)
vk0=var(noinfprobit0)
simk0=rmvnorm(10^5,mk0,2*vk0)
usk0=probitnoinflpost(simk0,y,X[,3:4])-dmvnorm(simk0,mk0,2*vk0,log=T)
bf0probit=mean(exp(usk))/mean(exp(usk0))
rm(mk0,vk0,simk0,usk0)
rm(noinfprobit0)

noinfprobit1=hmnoinfprobit(10000,y,X[,2:4],1)
mk1=apply(noinfprobit1,2,mean)
vk1=var(noinfprobit1)
simk1=rmvnorm(10^5,mk1,2*vk1)
usk1=probitnoinflpost(simk1,y,X[,2:4])-dmvnorm(simk1,mk1,2*vk1,log=T)
bf1probit=mean(exp(usk))/mean(exp(usk1))
rm(mk1,vk1,simk1,usk1)
rm(noinfprobit1)

noinfprobit2=hmnoinfprobit(10000,y,cbind(X[,1],X[,3:4]),1)
mk2=apply(noinfprobit2,2,mean)
vk2=var(noinfprobit2)
simk2=rmvnorm(10^5,mk2,2*vk2)
usk2=probitnoinflpost(simk2,y,cbind(X[,1],X[,3:4]))-dmvnorm(simk2,mk2,2*vk2,log=T)
bf2probit=mean(exp(usk))/mean(exp(usk2))
rm(mk2,vk2,simk2,usk2)
rm(noinfprobit2)

noinfprobit3=hmnoinfprobit(10000,y,cbind(X[,1:2],X[,4]),1)
mk3=apply(noinfprobit3,2,mean)
vk3=var(noinfprobit3)
simk3=rmvnorm(10^5,mk3,2*vk3)
usk3=probitnoinflpost(simk3,y,cbind(X[,1:2],X[,4]))-dmvnorm(simk3,mk3,2*vk3,log=T)
bf3probit=mean(exp(usk))/mean(exp(usk3))
rm(mk3,vk3,simk3,usk3)
rm(noinfprobit3)

noinfprobit4=hmnoinfprobit(10000,y,X[,1:3],1)
mk4=apply(noinfprobit4,2,mean)
vk4=var(noinfprobit4)
simk4=rmvnorm(10^5,mk4,2*vk4)
usk4=probitnoinflpost(simk4,y,X[,1:3])-dmvnorm(simk4,mk4,2*vk4,log=T)
bf4probit=mean(exp(usk))/mean(exp(usk4))
rm(mk4,vk4,simk4,usk4)
rm(noinfprobit4)
rm(simk,usk)

log10(bf0probit)
mkprob
diag(vkprob)
log10(c(bf1probit,bf2probit,bf3probit,bf4probit))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 4.4 Logit model

flatlogit=hmflatlogit(10000,y,X,1)
par(mfrow=c(4,3),mar=1+c(1.5,1.5,1.5,1.5))
plot(flatlogit[,1],type="l",xlab="Iterations",ylab=expression(beta[1]))
hist(flatlogit[1001:10000,1],nclass=50,prob=T,main="",xlab=expression(beta[1]))
acf(flatlogit[1001:10000,1],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(flatlogit[,2],type="l",xlab="Iterations",ylab=expression(beta[2]))
hist(flatlogit[1001:10000,2],nclass=50,prob=T,main="",xlab=expression(beta[2]))
acf(flatlogit[1001:10000,2],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(flatlogit[,3],type="l",xlab="Iterations",ylab=expression(beta[3]))
hist(flatlogit[1001:10000,3],nclass=50,prob=T,main="",xlab=expression(beta[3]))
acf(flatlogit[1001:10000,3],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(flatlogit[,4],type="l",xlab="Iterations",ylab=expression(beta[4]))
hist(flatlogit[1001:10000,4],nclass=50,prob=T,main="",xlab=expression(beta[4]))
acf(flatlogit[1001:10000,4],lag=1000,main="",ylab="Autocorrelation",ci=F)

S=readline(prompt="Type  <Return>   to continue : ")

mean(flatlogit[1001:10000,1])
mean(flatlogit[1001:10000,2])
mean(flatlogit[1001:10000,3])
mean(flatlogit[1001:10000,4])

exp(-2.5888*214.9+1.9967*130.1+2.1260*129.9+2.1879*9.5)/(1+exp(-2.5888*214.9+
1.9967*130.1+2.1260*129.9+2.1879*9.5))

S=readline(prompt="Type  <Return>   to continue : ")

noinflogit=hmnoinflogit(10000,y,X,1)
par(mfrow=c(4,3),mar=1+c(1.5,1.5,1.5,1.5))
plot(noinflogit[,1],type="l",xlab="Iterations",ylab=expression(beta[1]))
hist(noinflogit[1001:10000,1],nclass=50,prob=T,main="",xlab=expression(beta[1]))
acf(noinflogit[1001:10000,1],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(noinflogit[,2],type="l",xlab="Iterations",ylab=expression(beta[2]))
hist(noinflogit[1001:10000,2],nclass=50,prob=T,main="",xlab=expression(beta[2]))
acf(noinflogit[1001:10000,2],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(noinflogit[,3],type="l",xlab="Iterations",ylab=expression(beta[3]))
hist(noinflogit[1001:10000,3],nclass=50,prob=T,main="",xlab=expression(beta[3]))
acf(noinflogit[1001:10000,3],lag=1000,main="",ylab="Autocorrelation",ci=F)
plot(noinflogit[,4],type="l",xlab="Iterations",ylab=expression(beta[4]))
hist(noinflogit[1001:10000,4],nclass=50,prob=T,main="",xlab=expression(beta[4]))
acf(noinflogit[1001:10000,4],lag=1000,main="",ylab="Autocorrelation",ci=F)

mean(noinflogit[1001:10000,1])
mean(noinflogit[1001:10000,2])
mean(noinflogit[1001:10000,3])
mean(noinflogit[1001:10000,4])

S=readline(prompt="Type  <Return>   to continue : ")

mklog=apply(noinflogit,2,mean)
vklog=var(noinflogit)
simk=rmvnorm(100000,mklog,2*vklog)
usk=logitnoinflpost(simk,y,X)-dmvnorm(simk,mklog,2*vklog,log=T)

noinflogit0=hmnoinflogit(10000,y,X[,3:4],1)
mk0=apply(noinflogit0,2,mean)
vk0=var(noinflogit0)
simk0=rmvnorm(100000,mk0,2*vk0)
usk0=logitnoinflpost(simk0,y,X[,3:4])-dmvnorm(simk0,mk0,2*vk0,log=T)
bf0logit=mean(exp(usk))/mean(exp(usk0))
rm(mk0,vk0,simk0,usk0)
rm(noinflogit0)

noinflogit1=hmnoinflogit(10000,y,X[,2:4],1)
mk1=apply(noinflogit1,2,mean)
vk1=var(noinflogit1)
simk1=rmvnorm(100000,mk1,2*vk1)
usk1=logitnoinflpost(simk1,y,X[,2:4])-dmvnorm(simk1,mk1,2*vk1,log=T)
bf1logit=mean(exp(usk))/mean(exp(usk1))
rm(mk1,vk1,simk1,usk1)
rm(noinflogit1)

noinflogit2=hmnoinflogit(10000,y,cbind(X[,1],X[,3:4]),1)
mk2=apply(noinflogit2,2,mean)
vk2=var(noinflogit2)
simk2=rmvnorm(100000,mk2,2*vk2)
usk2=logitnoinflpost(simk2,y,cbind(X[,1],X[,3:4]))-dmvnorm(simk2,mk2,2*vk2,log=T)
bf2logit=mean(exp(usk))/mean(exp(usk2))
rm(mk2,vk2,simk2,usk2)
rm(noinflogit2)

noinflogit3=hmnoinflogit(10000,y,cbind(X[,1:2],X[,4]),1)
mk3=apply(noinflogit3,2,mean)
vk3=var(noinflogit3)
simk3=rmvnorm(100000,mk3,2*vk3)
usk3=logitnoinflpost(simk3,y,cbind(X[,1:2],X[,4]))-dmvnorm(simk3,mk3,2*vk3,log=T)
bf3logit=mean(exp(usk))/mean(exp(usk3))
rm(mk3,vk3,simk3,usk3)
rm(noinflogit3)

noinflogit4=hmnoinfprobit(10000,y,X[,1:3],1)
mk4=apply(noinflogit4,2,mean)
vk4=var(noinflogit4)
simk4=rmvnorm(100000,mk4,2*vk4)
usk4=logitnoinflpost(simk4,y,X[,1:3])-dmvnorm(simk4,mk4,2*vk4,log=T)
bf4logit=mean(exp(usk))/mean(exp(usk4))
rm(mk4,vk4,simk4,usk4)
rm(noinflogit4)
rm(simk,usk)

log10(bf0logit)
mklog
diag(vklog)
log10(c(bf1logit,bf2logit,bf3logit,bf4logit))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 4.5

airqual=na.omit(airquality)
ozone=cut(airqual$Ozone,c(min(airqual$Ozone),median(airqual$Ozone),max(airqual$Ozone)),include.lowest=T)
month=as.factor(airqual$Month)
tempe=cut(airqual$Temp,c(min(airqual$Temp),median(airqual$Temp),max(airqual$Temp)),include.lowest=T)

counts=table(ozone,tempe,month)
is.array(counts)
ftable(ozone,tempe,month)

S=readline(prompt="Type  <Return>   to continue : ")

counts=as.vector(counts)
ozo=gl(2,1,20)
temp=gl(2,2,20)
mon=gl(5,4,20)
model1=glm(counts~ozo+temp+mon+ozo*temp+ozo*mon+temp*mon,family=poisson())
anova(model1)
summary(model1)

S=readline(prompt="Type  <Return>   to continue : ")

x1=rep(1,20)
lulu=rep(0,20)
x2=x3=x4=x5=x6=x7=x8=x9=lulu
x2[ozo==2]=1
x3[temp==2]=1
x4[mon==2]=1
x5[mon==3]=1
x6[mon==4]=1
x7[mon==5]=1
x8[ozo==2 & temp==2]=1
x9[ozo==2 & mon==2]=1
x10=x11=x12=x13=x14=x15=x16=lulu
x10[ozo==2 & mon==3]=1
x11[ozo==2 & mon==4]=1
x12[ozo==2 & mon==5]=1
x13[temp==2 & mon==2]=1
x14[temp==2 & mon==3]=1
x15[temp==2 & mon==4]=1
x16[temp==2 & mon==5]=1

X=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)
model2=glm(counts~-1+X,family=poisson())
summary(model2)

S=readline(prompt="Type  <Return>   to continue : ")

flatloglin=hmflatloglin(10000,counts,X,0.5)

dev.off()
par(mfrow=c(4,4),mar=1+c(1.5,1.5,1.5,1.5),cex=0.8)
for (i in 1:16) plot(flatloglin[,i],type="l",ylab="",xlab="Iterations")

S=readline(prompt="Type  <Return>   to continue : ")

par(mfrow=c(4,4),mar=1+c(1.5,1.5,1.5,1.5),cex=0.8)
for (i in 1:16) hist(flatloglin[1001:10000,i],nclass=30,main="",ylab="",xlab="")

S=readline(prompt="Type  <Return>   to continue : ")

par(mfrow=c(4,4),mar=1+c(1.5,1.5,1.5,1.5),cex=0.8)
for (i in 1:16) acf(flatloglin[1001:10000,i],main="")

S=readline(prompt="Type  <Return>   to continue : ")

apply(flatloglin[1001:10000,],2,mean)


S=readline(prompt="Type  <Return>   to continue : ")

noinfloglin=hmnoinfloglin(10000,counts,X,0.5)
par(mfrow=c(4,4),mar=1+c(1.5,1.5,1.5,1.5),cex=0.8)
for (i in 1:16) plot(noinfloglin[,i],type="l",ylab="",xlab="Iterations")

S=readline(prompt="Type  <Return>   to continue : ")

par(mfrow=c(4,4),mar=1+c(1.5,1.5,1.5,1.5),cex=0.8)
for (i in 1:16) hist(noinfloglin[1001:10000,i],nclass=30,main="",ylab="",xlab="")

S=readline(prompt="Type  <Return>   to continue : ")

par(mfrow=c(4,4),mar=1+c(1.5,1.5,1.5,1.5),cex=0.8)
for (i in 1:16) acf(noinfloglin[1001:10000,i],main="")

S=readline(prompt="Type  <Return>   to continue : ")

apply(noinfloglin[1001:10000,],2,mean)

S=readline(prompt="Type  <Return>   to continue : ")

mklog=apply(noinfloglin,2,mean)
vklog=var(noinfloglin)
simk=rmvnorm(100000,mklog,2*vklog)
usk=loglinnoinflpost(simk,counts,X)-dmvnorm(simk,mklog,2*vklog,log=T)

noinfloglin1=hmnoinfloglin(10000,counts,cbind(X[,1:7],X[,9:16]),0.5)
mk1=apply(noinfloglin1,2,mean)
vk1=var(noinfloglin1)
simk1=rmvnorm(100000,mk1,2*vk1)
usk1=loglinnoinflpost(simk1,counts,cbind(X[,1:7],X[,9:16]))-dmvnorm(simk1,mk1,2*vk1,log=T)
bf1loglin=mean(exp(usk))/mean(exp(usk1))
rm(mk1,vk1,simk1,usk1)
rm(noinfloglin1)

noinfloglin2=hmnoinfloglin(10000,counts,cbind(X[,1:8],X[,13:16]),0.5)
mk2=apply(noinfloglin2,2,mean)
vk2=var(noinfloglin2)
simk2=rmvnorm(100000,mk2,2*vk2)
usk2=loglinnoinflpost(simk2,counts,cbind(X[,1:8],X[,13:16]))-dmvnorm(simk2,mk2,2*vk2,log=T)
bf2loglin=mean(exp(usk))/mean(exp(usk2))
rm(mk2,vk2,simk2,usk2)
rm(noinfloglin2)

noinfloglin3=hmnoinfloglin(10000,counts,X[,1:12],0.5)
mk3=apply(noinfloglin3,2,mean)
vk3=var(noinfloglin3)
simk3=rmvnorm(100000,mk3,2*vk3)
usk3=loglinnoinflpost(simk3,counts,X[,1:12])-dmvnorm(simk3,mk3,2*vk3,log=T)
bf3loglin=mean(exp(usk))/mean(exp(usk3))
rm(mk3,vk3,simk3,usk3)
rm(noinfloglin3)

mklog
diag(vklog)
log10(c(bf1loglin,bf2loglin,bf3loglin))
