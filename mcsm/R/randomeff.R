randomeff=function(nsim=10^3,a=10,b=30){
#Gibbs sampler for Girl's and Boys energy intake#
#Here we take a1=a2=a3 and b1=b2=b3#
  x1=log(Energy$Girls)
  x2=log(Energy$Boys)
  n1=length(x1);n2=length(x2);n=n1+n2
  x1bar=mean(x1);x2bar=mean(x2);xbar=(n1*x1bar+n2*x2bar)/n
  k=2				#number of treatments
  a1=a2=a3=a;b1=b2=b3=b
#---------------Initial values---------------------------
  mu=rep(xbar, nsim)
  theta1=rep(x1bar,nsim);theta2=rep(x2bar,nsim)
  sigma2mu=sigma2=tau2=rep((var(x1)+var(x2))/2,nsim)
  mu0=xbar
#--------------------------------------------------------
  for (i in 2:nsim){
#---------------------theta--------------------------------
    B1=sigma2[i-1]/(sigma2[i-1]+n1*tau2[i-1])
    theta1[i]=rnorm(1,mean=B1*mu[i-1]+(1-B1)*x1bar, sd=sqrt(tau2[i-1]*B1))
    B2=sigma2[i-1]/(sigma2[i-1]+n2*tau2[i-1])
    theta2[i]=rnorm(1,mean=B2*mu[i-1]+(1-B2)*x2bar, sd=sqrt(tau2[i-1]*B2))
#----------------------mu----------------------------------
    B=tau2[i-1]/(k*sigma2mu[i-1]+tau2[i-1])
    m=B*mu0+(1-B)*(n1*theta1[i]+n2*theta2[i])/n
    mu[i]=rnorm(1, mean=m,sd=sqrt(sigma2mu[i-1]*B))
#---------------------sigma---------------------------------
    sh1=(n/2)+a1
    ra1=(1/2)*(sum((x1-theta1[i])^2)+sum((x2-theta2[i])^2))+b1
    sigma2[i]=1/rgamma(1,shape=sh1,rate=ra1)
#---------------------tau------------------------------------
    sh2=(k/2)+a2
    ra2=(1/2)*(sum((theta1[i]-mu[i])^2)+sum((theta2[i]-mu[i])^2))+b2
    tau2[i]=1/rgamma(1,shape=sh2,rate=ra2)
#---------------------sigma2mu------------------------------------
    sh3=(1/2)+a3
    ra3=(1/2)*(mu[i]-mu0)^2+b3
    sigma2mu[i]=1/rgamma(1,shape=sh3,rate=ra3)
  }
#----------------Plots---------------------------------
  par(mfrow=c(2,3),mar=c(4,4,1,1))
  top=max(mu);bot=min(mu)
  hist(mu[(nsim/10):nsim],col="grey",breaks=25,xlab="",main=expression(mu),fre=F)
  hist(theta1[(nsim/10):nsim],col="grey",breaks=25,xlab="",main=expression(theta[1]),fre=F)
  hist(theta2[(nsim/10):nsim],col="grey",breaks=25,xlab="",main=expression(theta[2]),fre=F)
  hist(sqrt(sigma2mu[(nsim/10):nsim]),col="sienna",breaks=50,xlim=c(.5,4),xlab="",main=expression(sigma[mu]),fre=F)
  hist(sqrt(tau2[(nsim/10):nsim]), col="sienna",breaks=25,xlim=c(.5,4),xlab="",main=expression(tau),fre=F)
  hist(sqrt(sigma2[(nsim/10):nsim]), col="sienna",breaks=25,xlim=c(.5,2),xlab="",main=expression(sigma),fre=F)
  }
