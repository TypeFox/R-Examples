hastings=function(nsim=10^3){
#Hastings random walk algorithm#

a=c(.1,1,10);na=length(a)
x=array(0,c(na,nsim))

for (i in 1:na){
  acc=0
  for (j in 2:nsim){
	y<-x[i,(j-1)]+runif(1,min=-a[i],max=a[i])
	r=min(exp(-.5*((y^2)-(x[i,(j-1)]^2))),1)
	u<-runif(1);acc=acc+(u<r)
	x[i,j]<-y*(u<r)+x[i,(j-1)]*(u>r)
	}
}
#---------------Plots-------------------------------
par(mfrow=c(3,na),mar=c(4,4,2,1))
for(i in 1:na) plot((nsim-500):nsim,x[i,(nsim-500):nsim],ty="l",lwd=2,xlab="Iterations",ylab="",
	main=paste("Rate",(length(unique(x[i,]))/nsim),sep=" "))
for(i in 1:na){
  hist(x[i,],freq=F,xlim=c(-4,4),ylim=c(0,.4),col='grey',ylab="",xlab="",breaks=35,main="")
  curve(dnorm(x),lwd=2,add=T)
  }
for(i in 1:na) acf(x[i,], main="")
}
