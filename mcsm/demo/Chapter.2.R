#Section 2.1.1 Uniform generation

nsim=10^4       #number of random numbers
x=runif(nsim)
x1=x[-nsim]     #vectors to plot
x2=x[-1]	#adjacent pairs
par(mfrow=c(1,3),mar=c(4,4,1,1))
hist(x,fre=FALSE,xlab="",ylab="",main="")
plot(x1,x2,xlab=expression(x[t]),ylab=expression(x[t+1]))
acf(x,ylab="")

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
set.seed(1)
print(runif(5))
set.seed(1)
print(runif(5))

S=readline(prompt="Type  <Return>   to continue : ")

#Section 2.1.2 Inverse transform

nsim=10^4                #number of random variables
U=runif(nsim) 
X=-log(U)                #transforms of uniforms
Y=rexp(nsim)             #exponentials from R
X11(h=3.5)
par(mfrow=c(1,2),mar=c(2,2,2,2))        #plots
hist(X,freq=F,main="Exp from Uniform",ylab="",xlab="",ncl=150,col="grey",xlim=c(0,8))
curve(dexp(x),add=T,col="sienna",lwd=2)
hist(Y,freq=F,main="Exp from R",ylab="",xlab="",ncl=150,col="grey",xlim=c(0,8))
curve(dexp(x),add=T,col="sienna",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()

#Section 2.2 General transforms

print(system.time(test1()));print(system.time(test2()))

S=readline(prompt="Type  <Return>   to continue : ")

#Section 2.2.2 Discrete distributions

print(system.time(test3()));print(system.time(test4()))

S=readline(prompt="Type  <Return>   to continue : ")

#Section 2.2.3 Mixture representation

Nsim=10^4;n=6;p=.3
y=rgamma(nsim,n,rate=p/(1-p))
x=rpois(nsim,y)
hist(x,main="",freq=F,col="grey",breaks=40,xlab="",ylab="")
lines(1:50,dnbinom(1:50,n,p),lwd=2,col="sienna")

S=readline(prompt="Type  <Return>   to continue : ")

#Section 2.3 Mixture representation

print(optimize(f=function(x){dbeta(x,2.7,6.3)},interval=c(0,1),max=T)$objective)
print(optimize(f=function(x){dbeta(x,2.7,6.3)/dbeta(x,2,6)},max=T,interval=c(0,1))$objective)

betagen(Nsim=2500) #function provided in mcsm

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()
