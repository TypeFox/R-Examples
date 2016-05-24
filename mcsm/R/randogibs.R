randogibs=function(T=10^3){
#Random effect probit model from McCullach (1997)
beta0=-3
sigma0=1

#number of individuals
n=20
#number of replicas
m=35
#simulated data
x=matrix(sample(c(-1,0,1),n*m,rep=T),nrow=n)
y=matrix(0,ncol=m,nrow=n)
tru=rnorm(n)
for (i in 1:n)
  y[i,]=(runif(m)<1/(1+exp(-beta0*x[i,]-tru[i])))

#reference value w/o random effects
mlan=as.numeric(glm(as.vector(y)~as.vector(x)-1,family=binomial)$coef)

likecomp=function(beta,sigma,u){

  sum(as.vector(y)*(beta*as.vector(x)+rep(u,m)-rep(tru,m)))-
     sum(log(1+exp(beta*as.vector(x)+rep(u,m))))-sum(u^2)/(2*sigma^2)-log(sigma)
  }

#Another function for MCMC
gu=function(mu,i,beta,sigma){

  sum((y[i,]*(beta*x[i,]+mu))-log(1+exp(beta*x[i,]+mu)))-0.5*mu^2/sigma^2
  }

pro=function(beta,u){

  exp(as.vector(y)*(beta*as.vector(x)+rep(u,m)))/(1+exp(as.vector(y)*(beta*as.vector(x)+rep(u,m))))^2
  }

#Simple MCMC run on this model
beta=mlan
sigma=factor=1
acpt=bcpt=0
u=rnorm(n)
samplu=matrix(u,nrow=n)

for (iter in 2:T){

   #generation of the random effects
   u=rnorm(n)
   for (i in 1:n){

      mu=samplu[i,iter-1]
      u[i]=factor*sigma[iter-1]*rnorm(1)+mu

      if (log(runif(1))>gu(u[i],i,beta[iter-1],sigma[iter-1])-gu(mu,i,beta[iter-1],sigma[iter-1])){
         u[i]=mu
         }else{acpt=acpt+1}
     }
   
   samplu=cbind(samplu,u)

   #generation of (beta.sigma)
   sigma=c(sigma,1/sqrt(2*rgamma(1,0.5*n)/sum(u^2)))

   tau=sigma[iter]/sqrt(sum(as.vector(x^2)*pro(beta[iter-1],u)))
   betaprop=beta[iter-1]+rnorm(1)*factor*tau
   if (log(runif(1))>likecomp(betaprop,sigma[iter],u)-likecomp(beta[iter-1],sigma[iter],u)){
        betaprop=beta[iter-1]
	}else{bcpt=bcpt+1}

   beta=c(beta,betaprop)

   if (iter>100){
     if (bcpt<.1*iter) factor=factor/3
     if (bcpt>.9*iter) factor=factor*3
     }
}

library(coda)
plot(mcmc(cbind(beta,sigma)))

browser()
cumuplot(mcmc(cbind(beta,sigma)))

list(beta=beta,sigma=sigma)
}
