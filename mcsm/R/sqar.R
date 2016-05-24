sqar=function(T=10^4,multies=FALSE,outsave=FALSE,npara=5){
#Artificial situation of an AR model observed with noise
library(coda)

#parameters
rho=.85
tau=.2

#ensures the posterior is bimodal and not centered near sqrt(yc)
xc=1
while (xc>0){
  xm=rnorm(1,mean=-2)
  xc=rho*xm+rnorm(1);xp=rho*xc+rnorm(1)
  yc=xc^2+tau*rnorm(1)
  }

#targets
ef=function(x){
    -.5*((xm*rho-x)^2+(x*rho-xp)^2+(yc-x^2)^2/tau^2 )}

eef=function(x){
     exp(ef(x))}

if (!multies){ 	#single chain

#Metropolis-Hastings
smpl=NULL

for (scale in c(.1,.8)){
   xmc=sqrt(abs(yc))
   for (t in 2:T){
     prop=xmc[t-1]+scale*rnorm(1)
     if (log(runif(1))>ef(prop)-ef(xmc[t-1])) 
	prop=xmc[t-1]

     xmc=c(xmc,prop)
     }
   smpl=cbind(smpl,xmc)
   }

#ks.test
M=10
par(mfrow=c(2,2),mar=c(4,3,1,1))
for (io in 1:2){

xmc=smpl[,io]
print(geweke.diag(mcmc(xmc)))
print(heidel.diag(mcmc(xmc)))

thn=jitter(xmc[seq(1,T,by=M)])
kst=NULL
for (m in seq(T/(10*M),T/M,le=100)) kst=c(kst,ks.test(thn[1:(m/2)],thn[(m/2)+(1:(m/2))])$p)

hist(xmc,pro=T,col="grey85",nclass=150,main="",ylab="",xlab="")
ordin=apply(as.matrix(seq(min(xmc),max(xmc),le=200)),1,eef)
lines(seq(min(xmc),max(xmc),le=200),
ordin*max(density(xmc)$y)/max(ordin),lwd=2,col="gold4")
plot(seq(1,T,le=100),kst,pch=19,cex=.5,xlab="Iterations",ylab="p-value")

if (io==1) S=readline(prompt="Type  <Return>   to continue : ")
}

if (outsave) smpl
}else{ 		#multiple chains

smpl=nrmlz=NULL
scale=.5

for (t in 1:npara){

   xmc=sample(c(-1,1),1)*sqrt(abs(yc))+rnorm(1)*scale
   for (t in 2:T){

     prop=xmc[t-1]+scale*rnorm(1)
     if (log(runif(1))>ef(prop)-ef(xmc[t-1]))
        prop=xmc[t-1]

     xmc=c(xmc,prop)
     }

   smpl=cbind(smpl,mcmc(xmc))
   }

par(mfrow=c(1,2),mar=c(4,2,1,1))
plot(smpl[,1],type="l",ylim=range(smpl),xlab="Iterations",ylab="",col=heat.colors(10)[5])
for (t in 2:5) lines(smpl[,t],col=heat.colors(10)[6-t])
plot(density(smpl,n=1024),main="",ylab="",xlab=paste("Bandwith",format(density(smpl)$b,dig=3),sep=" "),lwd=2)

#Illustration of mcmc.list
S=readline(prompt="Type  <Return>   to continue : ")
plot(mcmc.list(mcmc(smpl[,1]),mcmc(smpl[,2]),mcmc(smpl[,3]),mcmc(smpl[,4]),mcmc(smpl[,5])))

if (outsave) smpl
}
}
