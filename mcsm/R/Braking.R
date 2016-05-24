Braking=function(nsim=10^3){
#This does Metropolis-Hastings for the Braking data
data=cars
x=data[,1]
y=data[,2]
n=length(x)
x2=x^2
MSR=225
SSR=10809
bhat=coefficients(lm(y~x+x2))
#-----------------likelihood---------------------------------------------
loglike=function(a,b,c,s2){-(n/2)*log(s2)-sum((y-a-b*x-c*x2)^2)/(2*s2)}
#-----------------candidate---------------------------------------------
dcand=function(a,b,c,s2){dnorm(a,bhat[1],sd=15,log=TRUE)+dnorm(b,bhat[2],sd=2,log=TRUE)
		+dnorm(c,bhat[3],sd=.06,log=TRUE)-(n/2)*log(s2)-SSR/(2*s2)}
#-----------------MH algorithm--------------------------------------------
b1hat=array(bhat[1],dim=c(nsim,1));b2hat=array(bhat[2],dim=c(nsim,1))
b3hat=array(bhat[3],dim=c(nsim,1));s2hat=array(MSR,dim=c(nsim,1));
for (i in 2:nsim)  {
  bcand=c(rnorm(1,mean=bhat[1],sd=15),rnorm(1,mean=bhat[2],sd=2),rnorm(1,mean=bhat[3],sd=.06),
 	1/rgamma(1,n/2,rate=SSR/2))
  test=min(exp(loglike(bcand[1],bcand[2],bcand[3],bcand[4])-loglike(b1hat[i-1],b2hat[i-1],
     b3hat[i-1],s2hat[i-1])+dcand(b1hat[i-1],b2hat[i-1],b3hat[i-1],s2hat[i-1])-dcand(bcand[1],bcand[2],bcand[3],bcand[4])),1);
   rho<-(runif(1)<test)
   b1hat[i]=bcand[1]*rho+b1hat[i-1]*(1-rho);
   b2hat[i]=bcand[2]*rho+b2hat[i-1]*(1-rho);
   b3hat[i]=bcand[3]*rho+b3hat[i-1]*(1-rho);
   s2hat[i]=bcand[4]*rho+s2hat[i-1]*(1-rho);
}
#------------------Plot----------------------------
plot(x,b1hat[nsim]+b2hat[nsim]*x+b3hat[nsim]*x2,type="l",col="grey",xlab="",ylab="",ylim=c(0,120),lwd=2)
for (i in (nsim-500):nsim){
   lines(x,b1hat[i]+b2hat[i]*x+b3hat[i]*x2,col="grey",lwd=2)}
lines(x,bhat[1]+bhat[2]*x+bhat[3]*x2,col="red",lwd=2)
points(x,y,pch=19)
}
