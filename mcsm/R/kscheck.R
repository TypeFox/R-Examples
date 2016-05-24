kscheck=function(T=10^3,heidel=FALSE){
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

betasave=beta;lambdasave=lambda

if (!heidel){	#apply Kolmogorov-Smirnov test:

#single chain KS:
ks=NULL
M=10

for (t in seq(T/10,T,le=100)){
    beta1=beta[1:(t/2)]
    beta2=beta[(t/2)+(1:(t/2))]
    beta1=beta1[seq(1,t/2,by=M)]
    beta2=beta2[seq(1,t/2,by=M)]
    ks=c(ks,ks.test(beta1,beta2)$p)
    }

#dual chain KS:
oldbeta=beta[seq(1,T,by=M)]
olks=ks

lambda=jitter(data/time)
beta=rgamma(1,gamma+10*alpha)/(delta+sum(lambda))

for (t in 2:T){
  lambda=rbind(lambda,rgamma(10,data+alpha)/(time+beta[t-1]))
  beta=c(beta,rgamma(1,gamma+10*alpha)/(delta+sum(lambda[t])))
  }

beta=beta[seq(1,T,by=M)]
ks=NULL
for (t in seq((T/(10*M)),(T/M),le=100)) 
    ks=c(ks,ks.test(beta[1:t],oldbeta[1:t])$p)

par(mar=c(4,4,1,1),mfrow=c(2,1))
plot(seq(1,T,le=100),olks,pch=19,cex=.7,xlab="Iterations",ylab="p-value")
plot(seq(1,T,le=100),ks,pch=19,cex=.7,xlab="Iterations",ylab="p-value")

}else{ 	# apply Cramer-von Mises test:
heidel.diag(mcmc(beta))
geweke.diag(mcmc(beta))
}

list(beta=betasave,lambda=lambdasave)
}
