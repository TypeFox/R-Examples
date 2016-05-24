normbyde=function(nsim=10^3,a=3){
#This gives an MH algorithm normal from DE.
dl=function(x,a)(a/2)*exp(-a*abs(x))	#Laplace density
X1=X2=array(0,dim=c(nsim,1))
X1[1]=X2[1]=rnorm(1)		#initialize the chain
acc1=acc2=0
for (i in 2:nsim){
	Y1=ifelse(runif(1) > 0.5, 1, -1) * rexp(1);Y2=Y1/a	
	rho1=min(1, dnorm(Y1)*dl(X1[i-1],1)/(dnorm(X1[i-1])*dl(Y1,1)))
	rho2=min(1, dnorm(Y2)*dl(X2[i-1],a)/(dnorm(X2[i-1])*dl(Y2,a)))
	test1=(runif(1)<rho1);acc1=acc1+test1
	test2=(runif(1)<rho2);acc2=acc2+test2
	X1[i]=X1[i-1] + (Y1-X1[i-1])*test1
	X2[i]=X2[i-1] + (Y2-X2[i-1])*test2
}
print(c(acc1/nsim,acc2/nsim))

par(mfrow=c(1,3),mar=c(4,4,2,1))
plot(cumsum(X1)/(1:nsim),type="l",ylim=c(-.25,.25),lwd=2,xlab="Iterations",ylab="",col="blue")
lines(cumsum(X2)/(1:nsim),lwd=2,xlab="",ylab="")
acf(X1);acf(X2)

print(c(mean(X1[(nsim/2):nsim]),var(X1[(nsim/2):nsim])))
print(c(mean(X2[(nsim/2):nsim]),var(X2[(nsim/2):nsim])))
}
