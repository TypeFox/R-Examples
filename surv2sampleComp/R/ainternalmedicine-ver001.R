# rm(list=ls())

######################################################################

rmstreg=function(y, delta, x, arm, tau, type="difference")
{

if(type!="difference" && type!="ratio" && type!="lossratio")
  print("Type must be difference, ratio or lossratio.")

if(type=="difference" || type=="ratio" || type=="lossratio"){

n=length(y)
x=cbind(1, x)
p=length(x[1,])

y0=pmin(y, tau)
d0=delta
d0[y0==tau]=1

d10=d0[arm==1]
d00=d0[arm==0]
y10=y0[arm==1]
y00=y0[arm==0]
x1=x[arm==1,]
x0=x[arm==0,]
n1=length(d10)
n0=length(d00)


id1=order(y10)
y10=y10[id1]
d10=d10[id1]
x1=x1[id1,]

id0=order(y00)
y00=y00[id0]
d00=d00[id0]
x0=x0[id0,]

fitc1=survfit(Surv(y10, 1-d10)~1)
fitc0=survfit(Surv(y00, 1-d00)~1)

weights1=d10/rep(fitc1$surv, table(y10))
weights0=d00/rep(fitc0$surv, table(y00))

weights=c(weights1, weights0)

if(type=="difference")
  {fitt=lm(c(y10,y00)~rbind(x1, x0)-1, weights=weights)
   beta0=fitt$coef
   
   error1=y10-as.vector(x1%*%beta0)
   score1=x1*weights1*error1

   error0=y00-as.vector(x0%*%beta0)
   score0=x0*weights0*error0
   }

if(type=="ratio")
  {fitt=glm(c(y10,y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
   beta0=fitt$coef
   
   error1=y10-exp(as.vector(x1%*%beta0))
   score1=x1*weights1*error1

   error0=y00-exp(as.vector(x0%*%beta0))
   score0=x0*weights0*error0
   }


if(type=="lossratio")
  {fitt=glm(c(tau-y10,tau-y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
   beta0=fitt$coef
   
   error1=tau-y10-exp(as.vector(x1%*%beta0))
   score1=x1*weights1*error1

   error0=tau-y00-exp(as.vector(x0%*%beta0))
   score0=x0*weights0*error0
   }



kappa.arm1=matrix(0, n1, p)
for(i in 1:n1)
   {kappa1=score1[i,]

    kappa2=apply(score1[y10>=y10[i],,drop=F], 2, sum)*(1-d10[i])/sum(y10>=y10[i])

    kappa3=rep(0, p)
    for(k in 1:n1)
       { if(y10[k]<=y10[i])
           kappa3=kappa3+apply(score1[y10>=y10[k],,drop=F], 2, sum)*(1-d10[k])/(sum(y10>=y10[k]))^2
       }

    kappa.arm1[i,]=kappa1+kappa2-kappa3
    }


kappa.arm0=matrix(0, n0, p)
for(i in 1:n0)
   {kappa1=score0[i,]

    kappa2=apply(score0[y00>=y00[i],,drop=F], 2, sum)*(1-d00[i])/sum(y00>=y00[i])

    kappa3=rep(0, p)
    for(k in 1:n0)
       {if(y00[k]<=y00[i])
           kappa3=kappa3+apply(score0[y00>=y00[k],,drop=F], 2, sum)*(1-d00[k])/(sum(y00>=y00[k]))^2
       }

    kappa.arm0[i,]=kappa1+kappa2-kappa3
    }


if(type=="difference")
  {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
   A=t(x)%*%x
   varbeta=solve(A)%*%gamma%*%solve(A)
   }

if(type=="ratio" || type=="lossratio")
  {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
   A=t(x*exp(as.vector(x%*%beta0)))%*%x
   varbeta=solve(A)%*%gamma%*%solve(A)
   }  


if(type=="difference")
  {beta0=beta0
   se0=sqrt(diag(varbeta))
   z0=beta0/se0
   p0=1-pchisq(z0^2, 1)
   cilow=beta0-se0*1.96
   cihigh=beta0+se0*1.96
   result=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, "lower .95"=cilow, "upper .95"=cihigh)
   }

if(type=="ratio" || type=="lossratio")
  {beta0=beta0
   se0=sqrt(diag(varbeta))
   z0=beta0/se0
   p0=1-pchisq(z0^2, 1)
   r0=exp(beta0)
   cilow=exp(beta0-se0*1.96)
   cihigh=exp(beta0+se0*1.96)
   result=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, "exp(coef)"=exp(beta0), "lower .95"=cilow, "upper .95"=cihigh)
   }

if(p==2)
   rownames(result)=c("intercept", "x")

if(p>2)
   rownames(result)=c("intercept", colnames(x[,-1]))

return(result)
}
}


##################################################################################


rmstaug=function(y, delta, x, arm, tau, type="difference")
{
if(type!="difference" && type!="ratio" && type!="lossratio")
  print("Type must be difference, ratio or lossratio.")

if(type=="difference" || type=="ratio" || type=="lossratio"){

n=length(y)
x=as.matrix(x)
p=length(x[1,])
pi=mean(arm)


y0=pmin(y, tau)
d0=delta
d0[y0==tau]=1

d10=d0[arm==1]
d00=d0[arm==0]
y10=y0[arm==1]
y00=y0[arm==0]
x1=x[arm==1,,drop=F]
x0=x[arm==0,,drop=F]
n1=length(d10)
n0=length(d00)


id1=order(y10)
y10=y10[id1]
d10=d10[id1]
x1=x1[id1,,drop=F]

id0=order(y00)
y00=y00[id0]
d00=d00[id0]
x0=x0[id0,,drop=F]

fitc1=survfit(Surv(y10, 1-d10)~1)
fitc0=survfit(Surv(y00, 1-d00)~1)

weights1=d10/rep(fitc1$surv, table(y10))
weights0=d00/rep(fitc0$surv, table(y00))

weights=c(weights1, weights0)

if(type=="difference")
  {fitt=lm(c(y10,y00)~rep(c(1, 0), c(n1, n0)), weights=weights)
   beta0=fitt$coef
   
   error1=y10-beta0[1]-beta0[2]
   score1=cbind(1, rep(1, n1))*weights1*error1

   error0=y00-beta0[1]
   score0=cbind(1, rep(0, n0))*weights0*error0
   }

if(type=="ratio")
  {fitt=glm(c(y10,y00)~rep(c(1, 0), c(n1, n0)), family="quasipoisson", weights=weights)
   beta0=fitt$coef
   
   error1=y10-exp(beta0[1]+beta0[2])
   score1=cbind(1, rep(1, n1))*weights1*error1

   error0=y00-exp(beta0[1])
   score0=cbind(1, rep(0, n0))*weights0*error0
   }


if(type=="lossratio")
  {fitt=glm(c(tau-y10,tau-y00)~rep(c(1, 0), c(n1, n0)), family="quasipoisson", weights=weights)
   beta0=fitt$coef
   
   error1=tau-y10-exp(beta0[1]+beta0[2])
   score1=cbind(1, rep(1, n1))*weights1*error1

   error0=tau-y00-exp(beta0[1])
   score0=cbind(1, rep(0, n0))*weights0*error0
   }



kappa.arm1=matrix(0, n1, 2)
for(i in 1:n1)
   {kappa1=score1[i,]

    kappa2=apply(score1[y10>=y10[i],,drop=F], 2, sum)*(1-d10[i])/sum(y10>=y10[i])

    kappa3=rep(0, 2)
    for(k in 1:n1)
       { if(y10[k]<=y10[i])
           kappa3=kappa3+apply(score1[y10>=y10[k],,drop=F], 2, sum)*(1-d10[k])/(sum(y10>=y10[k]))^2
       }

    kappa.arm1[i,]=kappa1+kappa2-kappa3
    }


kappa.arm0=matrix(0, n0, 2)
for(i in 1:n0)
   {kappa1=score0[i,]

    kappa2=apply(score0[y00>=y00[i],,drop=F], 2, sum)*(1-d00[i])/sum(y00>=y00[i])

    kappa3=rep(0, 2)
    for(k in 1:n0)
       {if(y00[k]<=y00[i])
           kappa3=kappa3+apply(score0[y00>=y00[k],,drop=F], 2, sum)*(1-d00[k])/(sum(y00>=y00[k]))^2
       }

    kappa.arm0[i,]=kappa1+kappa2-kappa3
    }


if(type=="difference")
  {
   A=cbind(c(n1+n0, n1), c(n1, n1))
   betainf=solve(A)%*%t(rbind(kappa.arm1, kappa.arm0))
   }

if(type=="ratio" || type=="lossratio")
  {mu1=exp(beta0[1]+beta0[2])
   mu0=exp(beta0[1])
   A=cbind(c(n1*mu1+n0*mu0, n1*mu1), c(n1*mu1, n1*mu1))
   
betainf=solve(A)%*%t(rbind(kappa.arm1, kappa.arm0))
   }  


aug=rbind(x1*(1-pi), -x0*pi)

fit=lm(t(betainf)~aug-1)

if(type=="difference")
  {beta0=beta0
   se0=sqrt(diag(betainf%*%t(betainf)))
   z0=beta0/se0
   p0=1-pchisq(z0^2, 1)
   cilow=beta0-se0*1.96
   cihigh=beta0+se0*1.96
   result.ini=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, "lower .95"=cilow, "upper .95"=cihigh)

   beta.aug=beta0-apply(aug%*%fit$coef,2,sum)
   se.aug=sqrt(diag(t(fit$res)%*%fit$res))*sqrt((n1+n0)/(n1+n0-p))
   z.aug=beta.aug/se.aug
   p.aug=1-pchisq(z.aug^2, 1)
   cilow.aug=beta.aug-se.aug*1.96
   cihigh.aug=beta.aug+se.aug*1.96
   result.aug=cbind(coef=beta.aug, "se(coef)"=se.aug, z=z.aug, p=p.aug, "lower .95"=cilow.aug, "upper .95"=cihigh.aug)
   }

if(type=="ratio" ||  type=="lossratio")
  {beta0=beta0
   se0=sqrt(diag(betainf%*%t(betainf)))
   z0=beta0/se0
   p0=1-pchisq(z0^2, 1)
   cilow=beta0-se0*1.96
   cihigh=beta0+se0*1.96
   result.ini=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, "exp(coef)"=exp(beta0), "lower .95"=exp(cilow), "upper .95"=exp(cihigh))

   beta.aug=beta0-apply(aug%*%fit$coef,2,sum)
   se.aug=sqrt(diag(t(fit$res)%*%fit$res))*sqrt((n1+n0)/(n1+n0-p))
   z.aug=beta.aug/se.aug
   p.aug=1-pchisq(z.aug^2, 1)
   cilow.aug=beta.aug-se.aug*1.96
   cihigh.aug=beta.aug+se.aug*1.96
   result.aug=cbind(coef=beta.aug, "se(coef)"=se.aug, z=z.aug, p=p.aug, "exp(coef)"=exp(beta.aug), "lower .95"=exp(cilow.aug), "upper .95"=exp(cihigh.aug))

   }

rownames(result.ini)=rownames(result.aug)=c("intercept", "arm")

return(list(result.ini=result.ini, result.aug=result.aug))
}
}


##################################################################################

# # 
# data=read.table("c:/temp/data-e4a03-covs.csv", sep=",", head=T, na.string="NA")
# data=data[,-1]
# data$stage[is.na(data$stage)]=2
# y=data$time
# delta=data$status

# arm=data[,3]
# tau0=40

# x=as.matrix(data[,c(3, 4,7,8)])
# rmstreg(y, delta, x, arm, tau=tau0, type="difference")
# rmstreg(y, delta, x, arm, tau=tau0, type="ratio")
# rmstreg(y, delta, x, arm, tau=tau0, type="lossratio")

# x=as.matrix(data[, c(4,7,8)])
# rmstaug(y, delta, x, arm, tau=tau0, type="difference")
# rmstaug(y, delta, x, arm, tau=tau0, type="ratio")
# rmstaug(y, delta, x, arm, tau=tau0, type="lossratio")


# x=arm
# rmstreg(y, delta, x, arm, tau=tau0, type="difference")
# rmstreg(y, delta, x, arm, tau=tau0, type="ratio")
# rmstreg(y, delta, x, arm, tau=tau0, type="lossratio")

# x=data[,4]
# rmstaug(y, delta, x, arm, tau=tau0, type="difference")
# rmstaug(y, delta, x, arm, tau=tau0, type="ratio")
# rmstaug(y, delta, x, arm, tau=tau0, type="lossratio")



######################################################################################
######################################################################################



