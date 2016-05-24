pimamh=function(Niter=10^4,scale=.01){
#Langevin implementation of the MCMC approximation to the marginal posterior
library(MASS) #to get Pima.tr
da=cbind(Pima.tr$type,Pima.tr$bmi)
da[,1]=da[,1]-1

like=function(a,b){
 sum(pnorm(q=a+outer(X=b,Y=da[,2],FUN="*"),log=T)*da[,1]+pnorm(q=-a-outer(X=b,Y=da[,2],FUN="*"),log=T)*(1-da[,1]))}

grad=function(a,b){
  don=pnorm(q=a+outer(X=b,Y=da[,2],FUN="*"))
  x1=sum((dnorm(x=a+outer(X=b,Y=da[,2],FUN="*"))/don)*da[,1]-
	(dnorm(x=-a-outer(X=b,Y=da[,2],FUN="*"))/(1-don))*(1-da[,1]))
  x2=sum((dnorm(x=a+outer(X=b,Y=da[,2],FUN="*"))/don)*da[,2]*da[,1]-
	(dnorm(x=-a-outer(X=b,Y=da[,2],FUN="*"))/(1-don))*da[,2]*(1-da[,1]))
  return(c(x1,x2))
  }

#starts from mle
the=matrix(glm(da[,1]~da[,2],family=binomial(link="probit"))$coef,ncol=2)
curmean=the[1,]+0.5*scale^2*grad(the[1,1],the[1,2])
likecur=like(the[1,1],the[1,2])

for (t in 2:Niter){

  prop=curmean+scale*rnorm(2)
  propmean=prop+0.5*scale^2*grad(prop[1],prop[2])
 
  if (log(runif(1))>like(prop[1],prop[2])-likecur-sum(dnorm(prop,mean=curmean,sd=scale,lo=T))+
     sum(dnorm(the[t-1,],mean=propmean,sd=scale,lo=T))){

     prop=the[t-1,];propmean=curmean
     }

   the=rbind(the,prop)
   curmean=propmean
}

# checking for the range
be1=seq(min(the[,1]),max(the[,1]),le=100)
be2=seq(min(the[,2]),max(the[,2]),le=130)
li=matrix(0,ncol=130,nro=100)
for (i in 1:100) for (j in 1:130) li[i,j]=like(be1[i],be2[j])

par(mar=c(4,4,1,1))
image(be1,be2,li,xlab=expression(beta[1]),ylab=expression(beta[2]))
contour(be1,be2,li,add=T,ncla=100)
subs=seq(1,Niter,le=10^3)
points(unique(the[subs,1]),unique(the[subs,2]),cex=.4,pch=19)
}
