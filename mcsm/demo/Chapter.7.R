#Section 7.2, Beta-binomial
betabi=function(x,a,b,n){
   beta(x+a,n-x+b)/(beta(a,b)*beta(x+1,n-x+1)*(n+1))}
nsim=10^4
n=15;a=3;b=7
X=T=rep(0,nsim)
T[1]=rbeta(1,a,b)               #initialize the chain
X[1]=rbinom(1,n,T[1])           #initialize the chain
for (i in 2:nsim){
        X[i]=rbinom(1,n,T[i-1])
        T[i]=rbeta(1,a+X[i],n-X[i]+b)
}
par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(X[2000:nsim],nclass=16,col="grey",freq=FALSE, xlim=c(0,15),main="",xlab="X")
curve(betabi(x,a,b,n),from=0, to=15,col="gold4",lwd=2,add=TRUE)
hist(T[2000:nsim],nclass=134,col="grey",freq=FALSE,xlim=c(0,.8),main="", xlab=expression(theta))
curve(dbeta(x,shape1=a,shape2=b),from=0, to=.8,col="sienna",lwd=2,add=TRUE)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.2, Energy example 7.3
data(Energy)
x=Energy$Girls
n=length(x);nsim=10^4
a=b=3;tau2=10;theta0=5
xbar=mean(x);sh1=(n/2)+a
sigma=theta=rep(0,nsim)                  #init arrays
sigma[1]=1/rgamma(1,shape=a,rate=b)      #init chains
B=sigma[1]/(sigma[1]+n*tau2)
theta[1]=rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
for (i in 2:nsim){
   B=sigma[i-1]/(sigma[i-1]+n*tau2)
   theta[i]=rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
   ra1=(1/2)*(sum((x-theta[i])^2))+b
   sigma[i]=1/rgamma(1,shape=sh1,rate=ra1)
   }
par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(log(theta[theta>0]),nclass=140,col="grey",freq=FALSE,main="",ylab="",
xlab=expression(theta),xlim=as.vector(quantile(log(theta[theta>0]),prob=c(.005,.995))))
hist(log(sigma),nclass=150,col="sienna",freq=FALSE,xlab=expression(sigma^2),
xlim=as.vector(quantile(log(sigma),prob=c(.005,.995))),main="",ylab="",)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.3, random effect example 7.5
randomeff()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.4, censored data model of Example 7.6
xdata=-c(3.64, 2.78, 2.91,2.85,2.54,2.62,3.16,2.21,4.05,2.19,2.97,4.32,
3.56,3.39,3.59,4.13,4.21,1.68,3.88,4.33)
m=length(xdata);n=30;a=3.5;nsim=10000
xbar=mean(xdata);that=rep(xbar,nsim)
zbar=rep(a,nsim)
for(i in 2:nsim){       
	temp=runif(n-m,min=pnorm(a,mean=that[i-1],sd=1),max=1)
        zbar[i]=mean(qnorm(temp,mean=that[i-1],sd=1))
        that[i]=rnorm(1,mean=(m/n)*xbar+(1-m/n)*zbar,sd=sqrt(1/n))}
par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(that[500:nsim],col="grey",breaks=25,,xlab=expression(theta),main="",freq=F)
hist(zbar[500:nsim],col="grey",breaks=25,main="",xlab= expression(bar(Z)),freq=F)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.4, multinomialmodel of Example 7.7
nsim=10^4
x=c(125,18,20,34);theta=rep(.5,nsim);z=rep(.5,nsim)
for (j in 2:nsim){
  theta[j]=rbeta(1,z[j-1]+x[4]+1,x[2]+x[3]+1)
  z[j]=rbinom(1,x[1],(theta[j]/(2+theta[j])))
  }
par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(theta[2000:nsim],freq=F,col="grey",breaks=25,main=expression(theta),xlab="")
hist(z[2000:nsim],freq=F,col="sienna",breaks=25,main="z",xlab="")

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.4, multinomial model of Example 7.8
nsim=5000
nA=186;nB=38;nAB=13;nO=284  		#Blood type data
pA=rep(.25,nsim);pB=array(.05,nsim)
for (i in 2:nsim){
  ZA=rbinom(1,nA,pA[i-1]^2/(pA[i-1]^2+2*pA[i-1]*(1-pA[i-1]-pB[i-1])))
  ZB=rbinom(1,nB,pB[i-1]^2/(pB[i-1]^2+2*pB[i-1]*(1-pA[i-1]-pB[i-1])))
  temp=rdirichlet(1,c(nA+nAB+ZA+1,nB+nAB+ZB+1,nA-ZA+nB-ZB+2*nO+1))
  pA[i]=temp[1];pB[i]=temp[2]
  }
par(mfrow=c(1,3),mar=c(4,4,2,1))
hist(pA,xlab=expression(p[A]),freq=F,col="grey",breaks=25,main="")
hist(pB,xlab=expression(p[B]),freq=F,col="sienna",breaks=25,main="")
hist(1-pA-pB,,xlab=expression(p[O]),freq=F,col="gold",breaks=25,main="")

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.4, mixture posterior of Example 7.9
gibbsmix()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.4, slice sampling in Example 7.10
fff=function(y)exp(-sqrt(y))/2;
nsim=10^3;x=rep(pi,nsim);x[1]=-log(runif(1))
for (i in 2:nsim){
  w=runif(1,min=0,max=fff(x[i-1]));x[i]=runif(1,min=0,max=(-log(2*w))^2)
}
X11(h=3.5);par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(x,xlab="Slice Sampler",main="",xlim=c(0,40),ylim=c(0,.25),freq=FALSE,col="grey",breaks=250)
curve(fff,add=T,lwd=3,col="sienna",xlim=c(0,40))
acf(x,xlab="Autocorrelation",ylab="",lwd=3)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.4, slice sampling in Example 7.11
data(challenger)
chares=challenge()
plot(chares$a,chares$b,type="l",xlab="a",ylab="b")

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.5, pump failure model  in Example 7.12
xdata=c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
Time=c(94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.05, 1.05, 2.10, 10.48)
nx=length(xdata)
nsim=10^4;alpha = 1.8;gamma= 0.01;delta=1
lambda=array(xdata*Time/sum(Time),dim=c(nsim,nx))
beta=rep(gamma*delta,nsim)
for (i in 2:nsim){
   for (j in 1:nx){
      lambda[i,j]=rgamma(1,shape=xdata[j]+alpha,rate=Time[j]+beta[i-1])
      }
   beta[i]=rgamma(1,shape=gamma+nx*alpha,rate=delta+sum(lambda[i,]))
   }
par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1.5));st=.5*nsim
hist(lambda[,1][st:nsim],breaks=25,col="grey",xlab="",main=expression(lambda[1]))
hist(lambda[,2][st:nsim],breaks=25,col="gold3",xlab="",main=expression(lambda[2]))
hist(beta[st:nsim],breaks=25,col="sienna",xlab="",main=expression(beta))
acf(lambda[,1][st:nsim],lwd=3,xlab="")
acf(lambda[,2][st:nsim],lwd=3,xlab="",col="gold3")
acf(beta[st:nsim],lwd=3,col="sienna",xlab="")

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.6, Reparameterization, Example 7.13
nsim=10^4
X=Y=rep(0,nsim)
rho=c(.3,.6,.9)
par(mfrow=c(1,3),mar=c(4,4,2,1))
for (t in 1:3){
  Std=sqrt(1-rho[t]^2)
  X[1]=rnorm(1)         #initialize the chain
  Y[1]=rnorm(1)         #initialize the chain
  for(i in 2:nsim){
     X[i]=rnorm(1,mean=rho[t]*Y[i-1],sd=Std)
     Y[i]=rnorm(1,mean=rho[t]*X[i],sd=Std)
     }
  acf(X,lwd=3)
  title(expression(rho==rho[t]),cex=2,col="sienna")
  }

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.6,1 Reparameterization, Example 7.14
reparareff()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.6.2, Rao-Blackwellization, Example 7.16
nsim=10^3
lam=RB=rep(313/360,nsim)
z=rep(0,13)
for (j in 2:nsim){
  top=round(lam[j -1]+6*sqrt(lam[j -1]))
  prob=dpois(c(4:top),lam[j -1])
  cprob=cumsum(prob/sum(prob))
  for(i in 1:13){z[i] = 4+sum(cprob<runif(1))}
  RB[j]=(313+sum(z))/360
  lam[j]=rgamma(1,360*RB[j],scale=1/360);
}
par(mfrow=c(1,3),mar=c(4,4,2,1))
hist(lam,col="grey",breaks=25,xlab="",main="Empirical Average")
plot(cumsum(lam)/c(1:nsim),ylim=c(1,1.05),type="l",lwd=1.5,ylab="")
lines(cumsum(RB)/c(1:nsim),col="sienna",lwd=1.5)
hist(RB,col="sienna",breaks=62,xlab="",main="Rao-Blackwell",xlim=c(1,1.05))

#Section 7.6,1 Reparameterization, Example 7.14
reparareff()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Section 7.6.3, improper priors
nsim=10^3
X=Y=rep(0,nsim)
X[1]=rexp(1)            #initialize the chain
Y[1]=rexp(1)            #initialize the chain
for (i in 2:nsim){
        X[i]=rexp(1,rate=Y[i-1])
        Y[i]=rexp(1,rate=X[i])
        }
st=0.1*nsim
par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(X,col="grey",breaks=25,xlab="",main="")
plot(cumsum(X)[(st+1):nsim]/c(1:(nsim-st)),type="l",lwd=1.5,ylab="")

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()
