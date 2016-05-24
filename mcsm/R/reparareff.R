reparareff=function(nsim=10^4,a=10,b=10){
#Gibbs sampler for Girl's and Boys energy intake#
#Compares Reparameterizations#
#Here we take a1=a2=a3 and b1=b2=b3#
x1=c(91,504,557,609,693,727,764,803,857,929,970,1043,1089,1195,1384,1713)
x2=c(457,645,736,790,899,991,1104,1154,1203,1320,1417,1569,1685,1843,2296,2710)
x1=log(x1)
x2=log(x2)
n1=length(x1);n2=length(x2);n=n1+n2
x1bar=mean(x1);x2bar=mean(x2);xbar=(n1*x1bar+n2*x2bar)/n
k=2				#number of treatments
a1=a2=a3=a;b1=b2=b3=b
#---------------Initial values---------------------------
nsim=10000
mu=array(xbar, c(nsim,1))
theta1=array(x1bar, c(nsim,1));theta2=array(x2bar, c(nsim,1))
sigma2mu=array((var(x1)+var(x2))/2,c(nsim,1))
sigma2=array((var(x1)+var(x2))/2,c(nsim,1))
tau2=array((var(x1)+var(x2))/2,c(nsim,1))
#---------------Reparameterization---------------------------
muR=array(xbar, c(nsim,1))
theta1R=array(x1bar, c(nsim,1));theta2R=array(x2bar, c(nsim,1))
sigma2muR=array((var(x1)+var(x2))/2,c(nsim,1))
sigma2R=array((var(x1)+var(x2))/2,c(nsim,1))
tau2R=array((var(x1)+var(x2))/2,c(nsim,1))
mu0=xbar
for(i in 2:nsim){
#---------------------theta--------------------------------
B1=sigma2[i-1]/(sigma2[i-1]+n1*tau2[i-1])
theta1[i]=rnorm(1,mean=B1*mu[i-1]+(1-B1)*x1bar, sd=sqrt(tau2[i-1]*B1))
B2=sigma2[i-1]/(sigma2[i-1]+n2*tau2[i-1])
theta2[i]=rnorm(1,mean=B2*mu[i-1]+(1-B2)*x2bar, sd=sqrt(tau2[i-1]*B2))
#---------------------repo--------------------------------
B1R=sigma2R[i-1]/(sigma2R[i-1]+n1*tau2R[i-1])
theta1R[i]=rnorm(1,mean=(1-B1R)*x1bar, sd=sqrt(tau2R[i-1]*B1R))
B2R=sigma2R[i-1]/(sigma2R[i-1]+n2*tau2R[i-1])
theta2R[i]=rnorm(1,mean=(1-B2R)*x2bar, sd=sqrt(tau2R[i-1]*B2R))
#----------------------mu----------------------------------
B=tau2[i-1]/(k*sigma2mu[i-1]+tau2[i-1])
m=B*mu0+(1-B)*(n1*theta1[i]+n2*theta2[i])/n
mu[i]=rnorm(1, mean=m,sd=sqrt(sigma2mu[i-1]*B))
#---------------------repo--------------------------------
BR=sigma2muR[i-1]/(sigma2muR[i-1]+sigma2R[i-1]/n)
m=(1-BR)*mu0+BR*(xbar-(n1*theta1R[i]+n2*theta2R[i])/n)
muR[i]=rnorm(1, mean=m,sd=sqrt(sigma2R[i-1]*BR/n))
#---------------------sigma---------------------------------
sh1=(n/2)+a1
ra1=(1/2)*(sum((x1-theta1[i])^2)+sum((x2-theta2[i])^2))+b1
sigma2[i]=1/rgamma(1,shape=sh1,rate=ra1)
#---------------------repo--------------------------------
ra1=(1/2)*(sum((x1-muR[i]-theta1R[i])^2)+sum((x2-muR[i]-theta2R[i])^2))+b1
sigma2[i]=1/rgamma(1,shape=sh1,rate=ra1)
#---------------------tau------------------------------------
sh2=(k/2)+a2
ra2=(1/2)*(sum((x1bar-theta1[i]-mu[i])^2)+sum((theta2[i]-mu[i])^2))+b2
tau2[i]=1/rgamma(1,shape=sh2,rate=ra2)
#---------------------repo--------------------------------
ra2=(1/2)*(theta1R[i]^2 + theta2R[i]^2)+b2
tau2[i]=1/rgamma(1,shape=sh2,rate=ra2)
#---------------------sigma2mu------------------------------------
sh3=(1/2)+a3
ra3=(1/2)*(mu[i]-mu0)^2+b3
sigma2mu[i]=1/rgamma(1,shape=sh3,rate=ra3)
}
#---------------------repo--------------------------------
ra3=(1/2)*(muR[i]-mu0)^2+b3
sigma2mu[i]=1/rgamma(1,shape=sh3,rate=ra3)
#----------------Plots---------------------------------
lim=c(-.05,.1)
par(mfrow=c(2,3),mar=c(4,4,4,1))
acf(mu,ylim=lim,lwd=2,main=expression(mu),cex=3)
acf(theta1,ylim=lim,lwd=2,main=expression(theta[1]))
acf(theta2,ylim=lim,lwd=2,main=expression(theta[2]))
acf(muR,ylim=lim,lwd=2,main=expression(mu))
acf(theta1R,ylim=lim,lwd=2,main=expression(theta[1]))
acf(theta2R,ylim=lim,lwd=2,main=expression(theta[2]))

print(cor(cbind(mu, theta1, theta2)))
print(cor(cbind(muR, theta1R, theta2R)))
}
