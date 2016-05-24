# Chapter 7 R commands

# Section 7.1

data(Eurostoxx50)
abnamro=Eurostoxx50[,2]
abnamro=ts(abnamro,freq=365-55*2,start=1998)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot.ts(abnamro,col="steelblue")
acf(abnamro,lag=365-55*2)
plot.ts(diff(abnamro),col="steelblue")
acf(diff(abnamro))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 7.2

p=10
T=260
dat=seqz=rnorm(T)
par(mfrow=c(2,2),mar=c(2,2,1,1))
for (i in 1:4){
  coef=runif(p,min=-.5,max=.5)
  for (t in ((p+1):T))
   seqz[t]=sum(coef*seqz[(t-p):(t-1)])+dat[t]
 plot(seqz,ty="l",col="sienna",lwd=2,ylab="")
 }

S=readline(prompt="Type  <Return>   to continue : ")

maxi=0
for (i in (1:10^4)) maxi=maxi+
    (max(Mod(polyroot(c(1,runif(10,-.5,.5)))))>1)
maxi/10^4

S=readline(prompt="Type  <Return>   to continue : ")

x=Eurostoxx50[, 4]
ar(x)
ar(x,method="ml")

S=readline(prompt="Type  <Return>   to continue : ")

niter=10^3
p=5
resAR5=ARmh(x=x,p=5,W=niter)
par(mfrow=c(3,3),mar=c(4,4,2,1))
hist(resAR5$ncomp,main="",xlab="p",ylab="",col="gold4",breaks=seq(-1,p+1,2))
par(new=TRUE)
plot(resAR5$ncomp,axes=F,cex=.3,xlab="",ylab="")
axis(side=4)
plot(resAR5$mus,type="l",col="steelblue4",xlab="Iterations",ylab=expression(mu))
plot(500:niter,resAR5$sigs[500:niter],type="l",col="steelblue4",xlab="Iterations",ylab=expression(sigma^2))
for (i in 1:min(3)) plot(500:niter,resAR5$psis[500:niter,i],type="l",col="steelblue4",xlab="Iterations",ylab=expression(psi))
plot(resAR5$llik,type="l",col="sienna4",xlab="Iterations",ylab="log-likelihood")

pst=matrix(1,ncol=6,nrow=niter)
pst[,1:5]=resAR5$psis
lame=apply(pst,1,polyroot)
plot((1/lame)[Mod(lame)>1],col="gold",cex=.3,xlab=expression(Re(lambda)),ylab=expression(Im(lambda)))
lines(seq(-1,1,.01),sqrt(1-seq(-1,1,.01)^2),col="sienna",lty=2,lwd=2)
lines(seq(-1,1,.01),-sqrt(1-seq(-1,1,.01)^2),col="sienna",lty=2,lwd=2)

pses=apply(resAR5$psis,2,mean)
mimine=mean(resAR5$mus)
predo=mimine-pses[1]*(x[p:(T-1)]-mimine)
for (i in 2:p) predo=predo-pses[i]*(x[(p-i+1):(T-i)]-mimine)
plot(x[(p+1):T],type="l",col="steelblue4",xlab="t",ylab="x")
lines(predo,lty=2,col="sienna4",lwd=1.8)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 7.3

x=Eurostoxx50[1:350,5]
arima(x,order=c(0,0,9))

S=readline(prompt="Type  <Return>   to continue : ")

niter=500
p=9
resMA9=MAmh(x=x,p=9,W=niter)

sube=function(x,n)
{
y=x
if (is.matrix(x))
{
t=dim(x)
if (t[1]>1000) y=y[seq(1,t[1],length=1000),]
if (t[2]>1000) y=y[,seq(1,t[2],length=1000)]
}
else y=y[seq(1,length(x),length=1000)]
y
}

start=100

par(mfrow=c(3,3),mar=c(4,4,2,1))
if (p>3)
{
hist(resMA9$ncomp,main="",xlab="p",ylab="",col="gold4",breaks=seq(-1,p+1,2));par(new=T);
plot(sube(1:niter),sube(resMA9$ncomp),axes=F,cex=.3,xlab="",ylab="")
axis(side=4)
}
if (p<4)
{
plot(sube(1:niter),sube(resMA9$ncomp),cex=.3,xlab="Iterations",ylab="Complex roots")
}
plot(sube(start:niter),sube(resMA9$mus[start:niter]),type="l",col="steelblue4",xlab="Iterations",ylab=expression(mu))
plot(sube(start:niter),sube(resMA9$sigs[start:niter]),type="l",col="steelblue4",xlab="Iterations",ylab=expression(sigma^2))
for (i in 1:min(3)) plot(sube(start:niter),sube(resMA9$psis[start:niter,i]),type="l",col="steelblue4",xlab="Iterations",ylab=expression(psi))
plot(sube(start:niter),sube(resMA9$llik[start:niter]),type="l",col="sienna4",xlab="Iterations",ylab="log-likelihood")

pst=matrix(1,ncol=(p+1),nrow=niter)
pst[,2:(p+1)]=-resMA9$psis
lame=apply(pst,1,polyroot)
plot(sube((1/lame)[Mod(lame)>1]),col="gold",cex=.3,xlab=expression(Re(lambda)),ylab=expression(Im(lambda)))
lines(seq(-1,1,.01),sqrt(1-seq(-1,1,.01)^2),col="sienna",lty=2,lwd=2)
lines(seq(-1,1,.01),-sqrt(1-seq(-1,1,.01)^2),col="sienna",lty=2,lwd=2)

plot(sube(start:niter),sube(resMA9$epsrec[start:niter,1]),ylim=range(resMA9$epsrec[start:niter,]),type="l",ylab=expression(epsilon),col="steelblue")
for (i in 2:p) lines(sube(start:niter),sube(resMA9$epsrec[start:niter,i]),col="steelblue",xlab="iterations")

#Section 7.5
data(Dnadataset)
res=hmhmm(M=10^2,Dnadataset$x)$par
par(mfrow=c(4,2),mar=c(3,3,2,1))
plot(res[,1],type="l",col="steelblue4",xlab=expression(p[11]))
plot(res[,2],type="l",col="steelblue4",xlab=expression(p[22]))
plot(res[,3],type="l",col="steelblue4",xlab=expression(q[1]^1))
plot(res[,4],type="l",col="steelblue4",xlab=expression(q[1]^2))
plot(res[,5],type="l",col="steelblue4",xlab=expression(q[2]^1))
plot(res[,6],type="l",col="steelblue4",xlab=expression(q[2]^2))
plot(res[,7],type="l",col="steelblue4",xlab=expression(q[3]^1))
plot(res[,8],type="l",col="steelblue4",xlab=expression(q[3]^2))
