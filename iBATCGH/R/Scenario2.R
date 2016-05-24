Scenario2 <-
function(sigmak=0.1){
g=100    #number of gene expression probes
m=1000   # number of CGH probes
s=100     # number of samples

#State specific mean and variance
mu1=-0.65
mu2=0
mu3=0.65
mu4=1.5

mu=c(mu1,mu2,mu3,mu4)

sigma1=0.1
sigma2=0.1
sigma3=0.1
sigma4=0.2

sigma=c(sigma1,sigma2,sigma3,sigma4)

#Generate xi
A=matrix(c(0.3,0.6,0.095,0.005,0.09,0.818,0.09,0.002,0.095,0.6,0.3,0.005,0.005,0.71,0.005,0.28),nrow=4,byrow=T)
AC=matrix(nrow=4,ncol=4)
for (i in 1:4){
AC[i,1]=A[i,1]
for (j in 2:4){
AC[i,j]=AC[i,(j-1)]+A[i,j]
}
}

A1=A
A1[1,]=c(0.7500, 0.1800, 0.0500, 0.020)
A1[2,]=c(0.4955, 0.0020, 0.4955, 0.007)
A1[3,]=c(0.0200, 0.1800, 0.7000, 0.100)
A1[4,]=c(0.0001, 0.3028, 0.1001, 0.597)

AC1=matrix(nrow=4,ncol=4)
for (i in 1:4){
AC1[i,1]=A1[i,1]
for (j in 2:4){
AC1[i,j]=AC1[i,(j-1)]+A1[i,j]
}
}

xi=matrix(2,nrow=s,ncol=m)
change=c(4:8,100:109,250:259,300,306,380,390,420:426,490:495,500:503,505,525:530,sample(531:1000,197))
change.complete=rep(0,m)
change.complete[change]=1
change.pos.two=which(change.complete==0)
change.partial=sample(change.pos.two[-1],375)
change.complete[change.partial]=2
q=10

for(j in 2:m){
if(change.complete[j]==1){
for(i in 1:s){
temp2=runif(1,0,1)
if(temp2<AC1[xi[i,j-1],1]){
xi[i,j]=1
}
if(AC1[xi[i,j-1],1]<=temp2 && temp2<AC1[xi[i,j-1],2]){
xi[i,j]=2
}
if(AC1[xi[i,j-1],2]<=temp2 && temp2<AC1[xi[i,j-1],3]){
xi[i,j]=3
}
if(AC1[xi[i,j-1],3]<=temp2){
xi[i,j]=4
}
}
}
if(change.complete[j]==2){
samples.to.change=sample(1:s,q)
for(i in 1:q){
temp2=runif(1,0,1)
if(temp2<AC1[xi[samples.to.change[i],j-1],1]){
xi[samples.to.change[i],j]=1
}
if(AC1[xi[samples.to.change[i],j-1],1]<=temp2 && temp2<AC1[xi[samples.to.change[i],j-1],2]){
xi[samples.to.change[i],j]=2
}
if(AC1[xi[samples.to.change[i],j-1],2]<=temp2 && temp2<AC1[xi[samples.to.change[i],j-1],3]){
xi[samples.to.change[i],j]=3
}
if(AC1[xi[samples.to.change[i],j-1],3]<=temp2){
xi[samples.to.change[i],j]=4
}

}
}

}

#Generate X

X=matrix(nrow=s,ncol=m)

for (i in 1:s){
  for(j in 1:m){
    X[i,j]=rnorm(1,mean=mu[xi[i,j]],sd=sigma[xi[i,j]])
  }
}

#Generate beta

beta=matrix(0,nrow=g,ncol=m)

beta[4,change[6:15]]=((-1)^(floor(runif(1,0,2))))*rnorm(10,mean=0.5,sd=0.3)
beta[10,change[16:25]]=((-1)^(floor(runif(1,0,2))))*rnorm(10,mean=0.5,sd=0.3)


#Generate epsilon
epsilon=NULL
for(i in 1:s){
epsilon=rbind(epsilon,rnorm(g,mean=0,sd=sigmak))
}
#Generate intercept
mu.g=rnorm(g,0,sd=0.1)
#Generate Y
Y=xi%*%t(beta)+mu.g+epsilon

##Empirical transition matrix

realA=Tran(xi)
realA[1,]=realA[1,]/sum(realA[1,])
realA[2,]=realA[2,]/sum(realA[2,])
realA[3,]=realA[3,]/sum(realA[3,])
realA[4,]=realA[4,]/sum(realA[4,])

##Beta different from zero
signbeta=which(beta!=0)

#Generate distances between probes
distance=rexp(m-1)
disfix=2*sum(distance)
return(list(Y=Y,X=X,Xi=xi,A=realA,mu=mu,Sd=sigma,coeff=beta,distance=distance,disfix=disfix))
}
