randogit=function(Tem=10^3,Tmc=10^2){
#Random effect logit model from McCullach (1997)
beta0=-3
sigma0=1

#number of individuals
n=20
#number of replicas
m=35
#simulated data
x=matrix(sample(c(-1,0,1),n*m,rep=TRUE),nrow=n)
y=matrix(0,ncol=m,nrow=n)
tru=rnorm(n)
for (i in 1:n)
  y[i,]=(runif(m)<1/(1+exp(-beta0*x[i,]-tru[i])))

#reference value w/o random effects
mlan=as.numeric(glm(as.vector(y)~as.vector(x)-1,family=binomial)$coef)

#MC replicas
omegas=exp(matrix(rnorm(Tmc*n),nrow=n))

#target function 
targ=function(beta){
   xs=exp(beta*x)
   xxs=x*xs
   losprod=0
   for (j in 1:m) for (t in 1:Tmc) 
     losprod=losprod+sum(xxs[,j]*omegas[,t]/(1+xs[,j]*omegas[,t]))

   losprod
   }

#dichotomous resolution of fixed point equation
dikoto=function(sol){

  prop=mlan
  val=targ(prop)
  scale=sigma0

  dir=1-2*(val>sol)#going up or down

  while (abs(val-sol)>10^(-3)){
  
     while (dir*(sol-val)>0){

         prop=prop+dir*scale
         val=targ(prop)
         }

     dir=1-2*(val>sol)
     scale=scale/10
     }

   prop
   }

#MC Likelihood function
likecomp=function(beta,sigma){

 exp(sum(as.vector(y)*(beta*as.vector(x)+rep(tru,m)))-sum(log(1+exp(beta*as.vector(x)+rep(tru,m)))))
}

#function for MCMC
gu=function(mu,i,beta,sigma){

  sum((y[i,]*(beta*x[i,]+mu))-log(1+exp(beta*x[i,]+mu)))-0.5*mu^2/sigma^2
  }

#MCEM iterations
beta=mlan
sigma=sigma0

diff=iter=factor=1
lval=likecomp(beta,sigma)

while (diff>10^(-2)){

  #simulating u's by MCMC
  samplu=matrix(tru,ncol=Tem,nrow=n)
  acpt=0

  for (t in 2:Tem){

    u=rnorm(n)

    for (i in 1:n){
    
      u[i]=factor*sigma[iter]*u[i]+samplu[i,t-1]

      if (log(runif(1))>gu(u[i],i,beta[iter],factor*sigma[iter])-gu(samplu[i,t-1],i,beta[iter],factor*sigma[iter])){
         u[i]=samplu[i,t-1]
         }else{acpt=acpt+1}

      }
      samplu[,t]=u
    }
   
   if (acpt<.1*Tem) factor=factor/3
   if (acpt>.9*Tem) factor=factor*3

  #EM equation
   sigma=c(sigma,sd(as.vector(samplu)))
   omegas=exp(samplu)
   beta=c(beta,dikoto(Tmc*sum(x*y)))
   lval=c(lval,likecomp(beta[iter+1],sigma[iter+1]))

   diff=max(abs(diff(beta[iter:(iter+1)])),abs(diff(sigma[iter:(iter+1)])))
   iter=iter+1
   Tem=Tem*2
   }

par(mfrow=c(2,1),mar=c(4,4,2,1))
plot(beta,sigma,pch=19,xlab=expression(beta),ylab=expression(sigma))
plot(lval,ty="l",xlab="Iterations",ylab=expression(L^c))
}
