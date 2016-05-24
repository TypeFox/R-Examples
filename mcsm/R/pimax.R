pimax=function(Nsim=10^3){
#Monte Carlo approximation to the marginal posterior
library(MASS) #to get Pima.tr
da=cbind(Pima.tr$type,Pima.tr$bmi)
da[,1]=da[,1]-1

like=function(a,b){
 sum(pnorm(q=a+outer(X=b,Y=da[,2],FUN="*"),log=T)*da[,1]+pnorm(q=-a-outer(X=b,Y=da[,2],FUN="*"),log=T)*(1-da[,1]))}

newlike=function(a,b){
 apply(pnorm(q=a+outer(X=b,Y=da[,2],FUN="*"),log=T)*da[,1],
 +pnorm(q=-a-outer(X=b,Y=da[,2],FUN="*"),log=T)*(1-da[,1]),1,sum)}

margap=function(a){
  b=rt(Nsim,df=5)
  dtb=dt(b,5,log=T)
  b=b*.1+.1
  themar=0
  for (i in 1:Nsim)
    themar=themar+exp(like(a,b[i])-dtb[i])
  }

#range of function values
aas=seq(-4,0,le=40)
raaj=matrix(0,ncol=40,nrow=100)
for (t in 1:100) raaj[t,]=apply(as.matrix(aas),1,margap) 

par(mfrow=c(3,1),mar=c(4,2,2,1))
plot(aas,raaj[sample(1:100,1),],type="l",col="white",lwd=2,ylim=c(0,max(raaj)),xlab=expression(theta[0]),ylab="")
polygon(c(aas,rev(aas)),c(apply(raaj,2,max),rev(apply(raaj,2,min))),col="gray80",border=F)
lines(aas,raaj[sample(1:100,1),],col="steelblue4",lwd=2)

thet=rt(Nsim,df=5);dtb=dt(thet,5,log=T);thet=thet*.1+.1

margbp=function(a){
  themar=0
  for (i in 1:Nsim)
    themar=themar+exp(like(a,thet[i])-dtb[i])
  }

rbbj=raaj*0
for (t in 1:100){

  thet=rt(Nsim,df=5);dtb=dt(thet,5,log=T);thet=thet*.1+.1
  rbbj[t,]=apply(as.matrix(aas),1,margbp)
  }

plot(aas,rbbj[sample(1:100,1),],type="l",col="white",lwd=2,ylim=c(0,max(rbbj)),xlab=expression(theta[0]),ylab="")
polygon(c(aas,rev(aas)),c(apply(rbbj,2,max),rev(apply(rbbj,2,min))),col="gray80",border=F)
lines(aas,rbbj[sample(1:100,1),],col="steelblue4",lwd=2)

plot(aas,apply(rbbj,2,mean),ty="l",col="sienna",lwd=2,xlab=expression(theta[0]),ylab="")
lines(aas,apply(raaj,2,mean),lty=2,col="gray30",lwd=2)
}
