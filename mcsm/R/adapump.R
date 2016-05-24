adapump=function(T=10^2,MM=10^2){
#Pump failure data
data=c(5,1,5,14,3,19,1,1,4,22)
time=c(94.32,15.72,62.88,125.76,5.24,31.44,1.05,1.05,2.10,10.48)
alpha=1.8;gamma=.01;delta=1

#Gibbs sampler
lambda=jitter(data/time)
beta=rgamma(1,sha=gamma+10*alpha)/(delta+sum(lambda))

for (t in 2:T){

  lambda=rbind(lambda,rgamma(10,sha=data+alpha)/(time+beta[t-1]))
  beta=c(beta,rgamma(1,sha=gamma+(10*alpha))/(delta+sum(lambda[t,])))
  }

cbeta=beta[T];clambda=lambda[T,]

post=function(lambda,beta){     # log posterior density
  sum(-(time+beta)*lambda+(data+alpha-1)*log(lambda))+
  dgamma(beta,shape=10*alpha+gamma,rate=delta,log=T)
  }

for (m in 1:MM){
  mu=c(apply(log(lambda),2,mean),mean(log(beta)))
  Sigma=.2*var(cbind(log(lambda),log(beta)))

  for (t in 1:T){
    prop=exp(rmunorm(mu,Sigma))

    if (log(runif(1))>post(prop[1:10],prop[11])-post(clambda,cbeta)+
    dmunorm(log(c(clambda,cbeta)),mu,Sigma,log=T)-dmunorm(log(prop),mu,Sigma,log=T)-
    sum(log(c(clambda,cbeta)))+sum(log(prop))) #jacobian
        prop=c(clambda,cbeta)

    clambda=prop[1:10]
    cbeta=prop[11]
    lambda=rbind(lambda,clambda)
    beta=c(beta,cbeta)
    }
}
burna=beta[(10^2+1):(10^2+MM*T)]
plot(burna,type="l",lwd=2,xlab="Iterations", ylab=expression(beta))

themus=thevar=NULL
for (m in 1:MM){

   themus=c(themus,rep(mean(beta[1:(T*m)]),T))
   thevar=c(thevar,rep(var(beta[1:(T*m)]),T))
   }
polygon(c(1:(MM*T),(MM*T):1),c(themus+2*sqrt(thevar),rev(themus-2*sqrt(thevar))),border=F,col="grey88")
lines(burna,lwd=2)
lines(themus,col="tomato",lwd=2)
}
