# Chapter 3 R commands

# Section 3.1

data(caterpillar)
y=log(caterpillar$y)
X=as.matrix(caterpillar[,1:8])
vnames=names(caterpillar)

par(mfrow=c(2,4),mar=c(4.2,2,2,1.2))
for (i in 1:8) plot(X[,i],y,xlab=vnames[i],pch=19,
col="sienna4",xaxt="n",yaxt="n")

S=readline(prompt="Type  <Return>   to continue : ")

# Section 3.2
X=scale(X)
summary(lm(y~X))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 3.3 Jeffreys' prior analysis

n=length(y)
p=dim(X)[2]
hatalpha=mean(y)
hatalpha

S=readline(prompt="Type  <Return>   to continue : ")

hatbeta=solve(t(X)%*%X)%*%t(X)%*%y
hatbeta

S=readline(prompt="Type  <Return>   to continue : ")

sigma2hat=sum((mean(y)+X%*%hatbeta-y)^2)/(33-9)
sum((mean(y)+X%*%hatbeta-y)^2)/(33-9-2)

S=readline(prompt="Type  <Return>   to continue : ")

varalpha=sigma2hat/n
varbeta=sigma2hat*diag(solve(t(X)%*%X))
HPD=matrix(0,9,2)
HPD[,1]=c(hatalpha,hatbeta)-qt(0.975,n-p-1)*sqrt(c(varalpha,varbeta))
HPD[,2]=c(hatalpha,hatbeta)+qt(0.975,n-p-1)*sqrt(c(varalpha,varbeta))
library(gplots)
dev.off()
par(mar=c(4,2,2,4))
plotCI(apply(HPD,1,mean),uiw=apply(HPD,1,diff)/2,pch=19,col="steelblue",lwd=2,xlab="index",yaxt="n")
axis(side=4,xlab=expression(hat(beta)))
abline(h=0,lwd=2,col="sienna",lty=2)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 2.4 Zellner

res1=BayesReg(y,X)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 3.6 Variable Selection

res2=ModChoBayesReg(y,X)
