#Section 6.3.1, Metropolis-Hasting definition, beta example 6.1

a=2.7; b=6.3; c=2.669 # initial values
nsim=5000
X=rep(runif(1),nsim)  # initialize the chain
for (i in 2:nsim){
   Y=runif(1)
   rho=dbeta(Y,a,b)/dbeta(X[i-1],a,b)
   X[i]=X[i-1] + (Y-X[i-1])*(runif(1)<rho)
   }
X11(h=3.5);plot(4500:4800,X[4500:4800],ty="l",lwd=2,xlab="Iterations",ylab="X")

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
ks.test(jitter(X),rbeta(5000,a,b))

par(mfrow=c(1,2),mar=c(2,2,1,1))
hist(X,nclass=150,col="grey",main="Metropolis-Hastings",fre=FALSE)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=TRUE)
hist(rbeta(5000,a,b),nclass=150,col="grey",main="Direct Generation",fre=FALSE)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=TRUE)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 6.2.2, comparison of gamma generators

a=4.85;nsim=10000;
X1=X2=array(0,dim=c(nsim,1))            #AR & MH
X1[1]=X2[1]=rgamma(1,a,rate=1)                #initialize the chain
for (i in 2:nsim){
    Y=rgamma(1,floor(a),rate=floor(a)/a)    #candidate
    rhoAR=(exp(1)*Y*exp(-Y/a)/a)^(a-floor(a))
    rhoMH=(dgamma(Y,a,rate=1)/dgamma(X2[i-1],a,rate=1))/(dgamma(Y,floor(a),
	rate=floor(a)/a)/dgamma(X2[i-1],floor(a),rate=floor(a)/a))
    rhoMH=min(rhoMH,1)
    X1[i]=Y*(runif(1)<rhoAR)                       #accepted values
    X2[i]=X2[i-1] + (Y-X2[i-1])*(runif(1)<rhoMH)
}
X1=X1[X1!=0]                   #The AR sample
par(mfrow=c(2,2),mar=c(4,4,2,2))
hist(X1,col="grey",nclas=125,freq=FALSE,xlab="",main="Accept-Reject",xlim=c(0,15))
curve(dgamma(x, a, rate=1),lwd=2,add=TRUE)
hist(X2[2500:nsim],nclas=125,col="grey",freq=FALSE,xlab="",main="Metropolis-Hastings",xlim=c(0,15))
curve(dgamma(x, a, rate=1),lwd=2,add=TRUE)
acf(X1,lag.max=50,lwd=2,col="red")              #Accept-Reject
acf(X2[2500:nsim],lag.max=50,lwd=2,col="blue")  #Metropolis-Hastings

mean(X1);var(X1);mean(X2[2500:nsim]);var(X2[2500:nsim])

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 6.2.2, Cauchy generator of example 6.2

nsim=10^4
X=Z=c(rt(1,1))  # initialize the chain from the stationary
for (t in 2:nsim){
    Y=rnorm(1) # candidate normal
    rho=min(1, dt(Y,1)*dnorm(X[t-1])/(dt(X[t-1],1)*dnorm(Y)))
    X[t]=X[t-1] + (Y-X[t-1])*(runif(1)<rho)
    Y=rt(1,.5)  # candidate t
    rho=min(1, dt(Y,1)*dt(Z[t-1],.5)/(dt(Z[t-1],1)*dt(Y,.5)))
    Z[t]=Z[t-1] + (Y-Z[t-1])*(runif(1)<rho)
    }

par(mfrow=c(3,2),mar=c(2,2,1,1))
plot(5000:5800,X[5000:5800],type="l",lwd=2,xlab="",ylab="")
plot(5000:5800,Z[5000:5800],type="l",lwd=2,xlab="",ylab="")
hist(X,nclass=100,col="grey",main="",xlab="",ylab="",fre=FALSE,xlim=c(-10,10))
curve(dt(x,1),col="sienna",lwd=2,add=TRUE)
hist(Z[abs(Z)<10],nclass=100,col="grey",main="",xlab="",ylab="",fre=FALSE,xlim=c(-10,10))
curve(dt(x,1),col="sienna",lwd=2,add=TRUE)
acf(X,lag.max=50,lwd=2,col="gold3");acf(Z,lag.max=50,lwd=2,col="gold3")

X11()
plot(cumsum(X<3)/(1:nsim),lwd=2,ty="l",ylim=c(.85,1),xlab="iterations",ylab="")
lines(cumsum(Z<3)/(1:nsim),lwd=2,col="sienna")

S=readline(prompt="Type  <Return>   to continue : ")
dev.off();dev.off()

#Section 6.2.2, Car braking set of example 6.3
Braking()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 6.2.3, Metropolis-Hastings on the original Hastings' example
hastings()

S=readline(prompt="Type  <Return>   to continue (warning: lengthy!): ")
dev.off()

#Section 6.2.3, Metropolis-Hastings on the mixture posterior
mhmix(Ni=10^4)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 6.2.4, Langevin Metropolis-Hastings on the probit posterior
pimamh()

S=readline(prompt="Type  <Return>   to continue (warning: lengthy!): ")
dev.off()

#Section 6.2.4, Langevin Metropolis-Hastings on the mixture posterior
mhmix(Ni=10^4,lange=TRUE,scale=.1)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 6.2.4, model selection
out=mochoice(10^4)
print(out$top)
apply(out$model,2,mean)

S=readline(prompt="Type  <Return>   to continue : ")

#Section 6.3, Acceptance raate
normbyde()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()
