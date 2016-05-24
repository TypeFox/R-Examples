#source("~/Work/AManal/MCMCglmm_2.19/vignettes/Figures/TEST.R")
library(MASS)
library(MCMCglmm)

verbose=FALSE
plotit=FALSE
DICtest=TRUE
SUMtest=TRUE
nsim<-10
nitt<-13000
thin<-10
burnin<-3000

psets<-c()

# poisson test 1.2 seconds OK
print("res1")
R<-diag(1)
res1<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1))
tpar<-c(1, 1)

psets<-c(psets, tpar)

for(i in 1:nsim){
l<-exp(rnorm(100,1,R))
y<-rpois(100,l)
data=data.frame(y1=y, y2=y)
if(DICtest){
m1<-MCMCglmm(cbind(y1)~1, family="poisson", data=data, prior=prior, verbose=verbose, nitt=3, thin=1, burnin=1, pl=T)
if(abs(-2*sum(dpois(data$y1, exp(m1$Liab[1,]), log=TRUE))-m1$Deviance[2])<1e-6){
 print("Deviance OK for univariate non-Gaussian (res1 - poisson)")
}else{
 stop("Deviance wrong for univariate non-Gaussian (res1 - poisson)")
}}

m1<-MCMCglmm(cbind(y1)~1, family="poisson", data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res1 different from expected"))
}
res1[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# multinomial test J=1 test 0.9 seconds  OK
print("res2")
res2<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1))
tpar<-c(1,1)
psets<-c(psets, tpar)
for(i in 1:nsim){
y<-rbinom(100,10,plogis(rnorm(100,1,1)))
data=data.frame(y1=y, y2=10-y)
m1<-MCMCglmm(cbind(y1,y2)~1, family="multinomial2", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res2 different from expected"))
}
res2[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# categorical test J=1 test 0.9 seconds OK
print("res3")
res3<-matrix(NA, nsim,3)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
tpar<-c(1,3/2, 1)
psets<-c(psets, tpar)
for(i in 1:nsim){
x<-rnorm(100)
y<-rbinom(100,1,plogis(rnorm(100,1+x*1.5,1)))
data=data.frame(y1=y, y2=1-y, x=x)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
m1<-MCMCglmm(y1~x, family="categorical", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res3 different from expected"))
}
res3[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# categorical test J=1 test with slice sampling 0.9 seconds OK
print("res3b")
res3b<-matrix(NA, nsim,2)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
tpar<-c(1,1)
psets<-c(psets, tpar)
for(i in 1:nsim){
y<-rbinom(100,1,plogis(rnorm(100,1,1)))
data=data.frame(y1=y, y2=1-y)
prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
if(DICtest){
m1<-MCMCglmm(y1~1, family="categorical", data=data, prior=prior,verbose=verbose, nitt=3, thin=1, burnin=1, pl=TRUE, slice=TRUE)
d<-sum(dbinom(data$y1, 1, plogis(m1$Liab[1,]), log=TRUE))
if(abs(-2*d-m1$Deviance[2])<1e-6){
 print("Deviance OK for sliced logit (res3b)")
}else{
 stop("Deviance wrong for sliced logit (res3b)")
}}
m1<-MCMCglmm(y1~1, family="categorical", data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin, slice=TRUE)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res3b different from expected"))
}
res3b[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# gauss with blocked random  1.5 seconds OK
print("res4")
res4<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior<-list(R=list(V=as.matrix(1), n=1), G=list(G1=list(V=G, n=1)))
tpar<-c(-1,1,2)
psets<-c(psets, tpar)
for(i in 1:nsim){
fac<-as.factor(sample(1:50,300,replace=TRUE))
ffac<-gl(2,150)
y<-mvrnorm(300, c(-1), R)+mvrnorm(50, c(0), G)[fac]
data=data.frame(y1=y, fac=fac, ffac=ffac)
if(DICtest){
m1<-MCMCglmm(y1~1,random=~fac, data=data, prior=prior, verbose=verbose, nitt=3, thin=1, burnin=1, pr=T)
if(abs(-2*sum(dnorm(data$y1, (cBind(m1$X, m1$Z)%*%m1$Sol[2,])@x, sqrt(m1$VCV[2,2]), log=TRUE))-m1$Deviance[2])<1e-6){
 print("Deviance OK for univariate Gaussian (res4)")
}else{
 stop("Deviance wrong for univariate Gaussian (res4)")
}}
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res4 different from expected"))
}
res4[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# binary with blocked random  1.5 seconds OK

print("res4c")
res4c<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
tpar<-c(-1,1,2)
psets<-c(psets, tpar)
prior=list(R=list(V=R, n=1, fix=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1), R)+mvrnorm(75, c(0), G)[fac]
data=data.frame(y1=rbinom(300, 1,plogis(y)), fac=fac, y2=y)
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior, family="categorical",verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res4c different from expected"))
}
res4c[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# gauss with correlated random 9.6 seconds OK

print("res5")
res5<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
tpar<-c(-1,1,2)
psets<-c(psets, tpar)
for(i in 1:nsim){
Ped<-cbind(1:400, c(rep(NA,100), sample(1:50,300,TRUE)),c(rep(NA,100), sample(51:100,300,TRUE)))
y<-mvrnorm(300, c(-1), R)+rbv(Ped,G)[101:400]
data=data.frame(y1=y, animal=as.factor(Ped[,1][101:400]), fac=gl(2,150))
system.time(m1<-MCMCglmm(y1~1, random=~animal, pedigree=Ped, data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin))
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res5 different from expected"))
}
res5[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# gauss with correlated random regression OK
print("res5b")
res5b<-matrix(NA, nsim,6)
R<-as.matrix(2)
G<-diag(2)
prior=list(R=list(V=R, nu=1),G=list(G1=list(V=G, nu=2)))
tpar<-c(-1, 1, 0,0,1,2)
psets<-c(psets, tpar)
for(i in 1:nsim){
Ped<-cbind(1:400, c(rep(NA,100), sample(1:50,300,TRUE)),c(rep(NA,100), sample(51:100,300,TRUE)))
x<-runif(300, -1,1)
bv<-rbv(Ped,G)[101:400,]
y<-mvrnorm(300, c(-1), R)+bv[,1]+bv[,2]*x
data=data.frame(y1=y, animal=as.factor(Ped[,1][101:400]), x=x)
m1<-MCMCglmm(y1~1, random=~us(1+x):animal, pedigree=Ped, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res5b different from expected"))
}
res5b[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# gauss with us random effect 3 seconds
print("res6")
res6<-matrix(NA, nsim,6)
G=matrix(c(1,0.5,0.5,2),2,2)
R<-as.matrix(2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
tpar<-c(-1, 1, 0.5, 0.5, 2, 2)
psets<-c(psets, tpar)
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
facf<-as.factor(sample(1:2,300,replace=TRUE))
rear<-sample(1:300)
fac<-gl(75,4)[rear]
facf<-gl(2,2,300)[rear]
y<-mvrnorm(300, c(-1), R)+rowSums(mvrnorm(75, c(0,0), G)[fac,]*cbind(facf==1,facf==2))
data=data.frame(y1=y, fac=fac, facf=facf)
m1<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res6 different from expected"))
}
res6[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# gauss with idh random effect 1.7 seconds
print("res7")
res7<-matrix(NA,nsim,4)
res7b<-matrix(NA,nsim,6)
R<-as.matrix(2)
G=matrix(c(1,0,0,2),2,2)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=2)))
tpar<-c(-1, 1, 2, 2)
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
facf<-as.factor(sample(1:2,300,replace=TRUE))
y<-mvrnorm(300, c(-1), R)+rowSums(mvrnorm(75, c(0,0), G)[fac,]*cbind(facf==1,facf==2))
data=data.frame(y1=y, fac=fac, facf=facf)
m1<-MCMCglmm(y1~1,random=~idh(facf):fac,  data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res7 different from expected"))
}
tpar<-c(-1, 1,0,0, 2, 2)
m2<-MCMCglmm(y1~1,random=~us(facf):fac,  data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m2$Sol, m2$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)){
print("res7b different from expected")
}
res7[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
res7b[i,]<-posterior.mode(mcmc(cbind(m2$Sol, m2$VCV)))
print(i)
}
psets<-c(psets, c(-1, 1, 2, 2), c(-1, 1,0,0, 2, 2))

# bivariate gauss with missing data 1.9 seconds
print("res8")
res8<-matrix(NA, nsim,6)
R=matrix(c(1,0.5,0.5,2),2,2)
prior=list(R=list(V=R, n=1))
tpar<-c(-1,1,1,0.5,0.5,2)
psets<-c(psets, tpar)
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=y[,2])
data$y1[sample(1:300, 150)]<-NA
if(DICtest){
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~us(trait):units, data=data, prior=prior, verbose=verbose, nitt=3, thin=1, burnin=1, pl=TRUE)
pred1<-(m1$X%*%m1$Sol[2,])[1:300]
pred2<-(m1$X%*%m1$Sol[2,])[301:600]
l1<-m1$Liab[1,1:300]
l2<-m1$Liab[1,301:600]
Ritt<-matrix(m1$VCV[2,],2,2)
cpred2<-pred2+(Ritt[1,2]/Ritt[1,1])*(l1-pred1)
cR<-Ritt[2,2]-(Ritt[1,2]^2)/Ritt[1,1]
d<-0
for(j in 1:300){
if(is.na(data$y1[j])){
d<-d+dnorm(data$y2[j], cpred2[j],sqrt(cR), log=TRUE)
}else{
d<-d+mvtnorm::dmvnorm(c(data$y1[j], data$y2[j]), c(pred1[j], pred2[j]), Ritt, log=TRUE) 
}
}
if(abs(-2*d-m1$Deviance[2])<1e-6){
 print("Deviance OK for bivariate Gaussian with missing data (res8)")
}else{
 stop("Deviance wrong for bivariate Gaussian with missing data (res8)")
}}
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin, pl=TRUE)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res8 different from expected"))
}
res8[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate gauss idh residual 1.7 seconds
print("res9")
res9<-matrix(NA, nsim,4)
R=matrix(c(1,0,0,2),2,2)
prior=list(R=list(V=R, n=2))
tpar<-c(-1,1,1,2)
psets<-c(psets, tpar)
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=y[,2])
data$y1[sample(1:300, 50)]<-NA
if(DICtest){
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior, verbose=verbose, nitt=3, thin=1, burnin=1)
d<-sum(dnorm(data$y1[which(!is.na(data$y1))], (m1$X%*%m1$Sol[2,])[1:300][which(!is.na(data$y1))], sqrt(m1$VCV[2,1]), log=TRUE))
d<-d+sum(dnorm(data$y2, (m1$X%*%m1$Sol[2,])[301:600], sqrt(m1$VCV[2,2]), log=TRUE))
if(abs(-2*d-m1$Deviance[2])<1e-6){
 print("Deviance OK for idh bivariate Gaussian with missing data (res9)")
}else{
 stop("Deviance wrong for idh bivariate Gaussian with missing data (res9)")
}}

m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res9 different from expected"))
}
res9[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate gauss + random 2.8 seconds
print("res10")
res10<-matrix(NA, nsim,6)
R=matrix(c(1,0,0,2),2,2)
G=matrix(c(2,0,0,1),2,2)
tpar<-c(-1,1,2,1,1,2)
psets<-c(psets, tpar)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1,1), R)+mvrnorm(75, c(0,0), G)[fac,]
data=data.frame(y1=y[,1], y2=y[,2], fac=fac)
m1<-MCMCglmm(cbind(y1,y2)~trait-1, random=~idh(trait):fac, family=c("gaussian","gaussian"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res10 different from expected"))
}
res10[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# bivariate binoimal + random 5.9 seconds
res11<-matrix(NA, nsim,8)
R=matrix(c(1,0,0,2),2,2)
G=matrix(c(2,0.5,0.5,1),2,2)
tpar<-c(-1,1,2,0.5,0.5,1,1,2)
psets<-c(psets, tpar)
prior=list(R=list(V=R, n=1),G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, c(-1,1), R)+mvrnorm(75, c(0,0), G)[fac,]
y1<-rbinom(300,10,plogis(y[,1]))
y2<-rbinom(300,10,plogis(y[,2]))
data=data.frame(y1s=y1,y1f=10-y1,y2s=y2,y2f=10-y2, fac=fac)
m1<-MCMCglmm(cbind(y1s,y1f,y2s,y2f)~trait-1, random=~us(trait):fac, family=c("multinomial2","multinomial2"), rcov=~idh(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res11 different from expected"))
}
res11[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}


# gauss with blocked random - all missing 1.9 seconds
print("res12")
res12<-matrix(NA, nsim,3)
R<-as.matrix(2)
G<-as.matrix(1)
prior=list(R=list(V=R, n=100),G=list(G1=list(V=G, n=100)), B=list(mu=0, V=0.000000001))
tpar<-c(0,1,2)
psets<-c(psets, tpar)
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
y<-rep(NA,300)
data=data.frame(y1=y, fac=fac)
m1<-MCMCglmm(y1~1,random=~fac,  data=data, prior=prior, singular.ok=TRUE,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res12 different from expected"))
}

res12[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

# Jennys data - once as gaussian once as binomial 32 seconds

if(file.exists("~/Work/Jenny/Data/Intermediate/ThirdC.R")){

  m1R<-dget("~/Work/Jenny/Data/Intermediate/ThirdC.R")
  res13<-matrix(NA, nsim,39)
  res14<-matrix(NA, nsim,39)
  print("res13")
  print("res14")

  firstP<-read.table("~/Work/Jenny/Data/Raw/Third_paternal.txt", header=T)
  firstP$day<-as.factor(firstP$day)

  coef<-apply(cbind(m1R$Sol, m1R$VCV), 2, median)
  G<-list(G1=matrix(coef[9:33],5,5), G2=as.matrix(coef[34]))


  R<-diag(coef[35:39])
  prior=list(R=list(V=R, n=5),G=list(G1=list(V=G[[1]], n=1), G2=list(V=G[[2]], n=1)))


  for(i in 1:nsim){
    location<-coef[firstP$virus]+c(0,coef[1])[((firstP$virus!=levels(firstP$virus)[1]))+1]+c(0,coef[6:8])[firstP$day]
    resid<-rnorm(dim(firstP)[1], 0, sqrt(diag(R)[firstP$virus]))
    reffects<-rowSums(mvrnorm(nlevels(firstP$line),rep(0,5),G[[1]])[as.numeric(firstP$line),]*outer(firstP$virus, levels(firstP$virus), "=="))
    reffects<-reffects+rnorm(nlevels(firstP$f2rep),0,sqrt(G[[2]]))[as.numeric(firstP$f2rep)]
    l<-location+resid+reffects
    firstP$l<-l
    y<-rbinom(dim(firstP)[1], firstP$total, plogis(l))
    firstP$noalive<-y
    firstP$nodead<-firstP$total-y
    tpar<-coef
    m1test<-MCMCglmm(l~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("gaussian"), data=firstP, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1test$Sol, m1test$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test$Sol, m1test$VCV)))[,2]<tpar)/length(tpar), "res13 different from expected"))
}


    m1test2<-MCMCglmm(cbind(noalive, nodead)~virus+day, random=~us(virus):line+f2rep, rcov=~idh(virus):units, family=c("multinomial2"), data=firstP, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m2)
}
if(plotit){
plot(mcmc(cbind(m1test2$Sol, m1test2$VCV)), ask=FALSE)
}

if(any(HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1test2$Sol, m1test2$VCV)))[,2]<tpar)/length(tpar), "res14 different from expected"))
}

    res13[i,]<-posterior.mode(mcmc(cbind(m1test$Sol, m1test$VCV)))
    res14[i,]<-posterior.mode(mcmc(cbind(m1test2$Sol, m1test2$VCV)))
    print(i)
  }
psets<-c(psets, coef, coef)

}else{
  print("file ~/Work/Jenny/Data/Intermediate/ThirdC.R does not exist: skipping res 13 and 14")
}

res15<-matrix(NA, nsim,3)
print("res15")
R<-as.matrix(2)
prior=list(R=list(V=R, n=1))
tpar<-c(-1,1,2)

for(i in 1:nsim){
mev<-rchisq(300,10)
y<-mvrnorm(300, c(-1), R)+rnorm(300, c(0), sqrt(mev))
data=data.frame(y1=y)
m1<-MCMCglmm(y1~1,data=data, prior=prior, mev=mev,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res15 different from expected"))
}
res15[i,]<-apply(mcmc(cbind(m1$Sol, m1$VCV)), 2, median)
print(i)
}
psets<-c(psets, tpar)


# censored gaussian data

res17a<-matrix(NA, nsim,2)
res17b<-matrix(NA, nsim,2)
print("res17")
for(i in 1:nsim){
l<-rnorm(500,0,1)
y<-cut(l,c(-Inf,seq(-4,4,1), Inf), , include.lowest=TRUE)
y1<-c(c(-Inf,seq(-4,4,1), Inf),Inf)[y]
y2<-c(c(-Inf,seq(-4,4,1), Inf),Inf)[as.numeric(y)+1]
ym<-(y1+y2)/2
tpar<-c(0,1)
data=data.frame(y1=y1, y2=y2, ym=ym)
if(any(y1==-Inf) | any(y2==Inf)){
  m1$VCV<-rep(-1,1000)
  m1$Sol<-rep(-1,1000)
}else{
  m1<-MCMCglmm(ym~1, data=data, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
}
  m2<-MCMCglmm(cbind(y1, y2)~1, data=data, family="cengaussian",verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m2)
}
  if(plotit){
    plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
  }
  if(plotit){
    plot(mcmc(cbind(m2$Sol, m2$VCV)), ask=FALSE)
  }
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res17a different from expected"))
}
if(any(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)/length(tpar), "res17b different from expected"))
}
res17a[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
res17b[i,]<-posterior.mode(mcmc(cbind(m2$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar, tpar)

# censored Possion data
res18<-matrix(NA, nsim,2)
#print("res18")
tpar<-c(0,1)
prior=list(R=list(V=diag(1), n=1))
for(i in 1:nsim){
l<-rpois(100, exp(rnorm(100,0,1)))
y<-cut(l,c(seq(0,10,2), Inf), include.lowest=TRUE, right=FALSE)
y1<-c(c(seq(0,10,2), Inf),Inf)[y]
y2<-c(c(seq(0,10,2), Inf),Inf)[as.numeric(y)+1]-1
y1[100]<-l[100]
y2[100]<-l[100]
data=data.frame(y1=y1, y2=y2)
m1<-MCMCglmm(cbind(y1, y2)~1, data=data, family="cenpoisson", prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res18 different from expected"))
}
res18[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
#print(i)
}

psets<-c(psets, tpar)

res19<-matrix(0, nsim, 6)
res19b<-matrix(0, nsim, 6)

print("res19")
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2))
tune=list(diag(2))
tpar<-c(1,-1,1,-1/2, 1, 1)
for(i in 1:nsim){

x<-rnorm(300)
l<-mvrnorm(300,c(1, -1), R)
l[,1]<-l[,1]+x*1
l[,2]<-l[,2]-x*1/2

y<-rbinom(300, 1, 1-plogis(l[,2]))
y[which(y==1)]<-rpois(sum(y==1), exp(l[,1][which(y==1)]))



data=data.frame(y1=y, x=x)
m1<-MCMCglmm(y1~trait+trait:x-1, rcov=~idh(trait):units, data=data, family="zipoisson",prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res19 different from expected"))
}
res19[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}

psets<-c(psets, tpar)

print("res19b")
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2),G=list(G1=list(V=1, n=1)))
tune=list(diag(2))
tpar<-c(1, -0.5, 0.2, 1,1,1)
    for(i in 1:nsim){
	x<-rnorm(300)
	l<-mvrnorm(300,c(0,-0.5), R)
        fac<-gl(50,6)
        r<-rnorm(50)
        y<-rbinom(300, 1, 1-plogis(l[,2]))
        y[which(y==1)]<-rpois(sum(y==1), exp(1+l[,1]+r[fac]+0.2*x)[which(y==1)])


	data=data.frame(y1=y, x=x, fac=fac)
	m1<-MCMCglmm(y1~trait+at.level(trait, 1):x-1, random=~us(at.level(trait,1)):fac, rcov=~idh(trait):units, data=data, family="zipoisson",prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
       if(plotit){
          plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
       }
       if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
          print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res19b different from expected"))
       }
       res19b[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
       print(i)
    }

   psets<-c(psets, tpar)


# bivariate gauss + categorical  residual 1.7 seconds

res20<-matrix(NA, nsim,6)
print("res20")
tune<-diag(2)
R=matrix(c(2,0.25,0.25,1),2,2)
prior=list(R=list(V=R, n=2, fix=2))
tpar<-c(-1,1,2,0.25, 0.25,1)
for(i in 1:nsim){
y<-mvrnorm(300, c(-1,1), R)
data=data.frame(y1=y[,1], y2=rbinom(300, 1, plogis(y[,2])), y3=y[,2])
if(DICtest){
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","categorical"), rcov=~us(trait):units, data=data, prior=prior, verbose=verbose, nitt=3, thin=1, burnin=1, pl=TRUE)
pred1<-(m1$X%*%m1$Sol[2,])[1:300]
pred2<-(m1$X%*%m1$Sol[2,])[301:600]
l1<-m1$Liab[1,1:300]
l2<-m1$Liab[1,301:600]
Ritt<-matrix(m1$VCV[2,],2,2)
cpred1<-pred1+(Ritt[1,2]/Ritt[2,2])*(l2-pred2)
cR<-Ritt[1,1]-(Ritt[1,2]^2)/Ritt[2,2]
d<-sum(dnorm(data$y1, cpred1, sqrt(cR), log=TRUE))
d<-d+sum(dbinom(data$y2, 1, plogis(l2), log=TRUE))
if(abs(-2*d-m1$Deviance[2])<1e-6){
 print("Deviance OK for bivariate Gaussian/categorical (res20)")
}else{
 stop("Deviance wrong for idh bivariate Gaussian/categorical (res20)")
}}
m1<-MCMCglmm(cbind(y1,y2)~trait-1, family=c("gaussian","categorical"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res20 different from expected"))
}
res20[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

# gauss random regression

res21<-matrix(NA, nsim,5)
res21b<-matrix(NA, nsim,5)
res21c<-matrix(NA, nsim,7)
print("res21")
G=matrix(c(2,0,0,1),2,2)
R=matrix(1,1,1)
prior=list(R=list(V=R, n=1), G=list(G1=list(V=G[1,1,drop=FALSE], n=1), G2=list(V=G[2,2,drop=FALSE], n=0.1)))
prior2=list(R=list(V=R, n=1), G=list(G1=list(V=G, n=1)))

int.slope<-mvrnorm(300, c(0,0), G)
ind<-gl(300,3)
time<-rnorm(900)
y<-int.slope[,1][ind]+time*int.slope[,2][ind]
y<-y+rnorm(900,-1,R)
data=data.frame(y1=y, time=time, ind=ind)

for(i in 1:nsim){
int.slope<-mvrnorm(300, c(0,0), G)
ind<-gl(300,3)
time<-rnorm(900)
y<-int.slope[,1][ind]+time*int.slope[,2][ind]
y<-y+rnorm(900,-1,R)
data=data.frame(y1=y, time=time, ind=ind)
	m1<-MCMCglmm(y1~time, random=~ind+us(time):ind, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
	m2<-MCMCglmm(y1~time, random=~idh(1+time):ind, data=data, prior=prior2,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
        m3<-MCMCglmm(y1~time, random=~us(1+poly(time,1, raw=TRUE)):ind, data=data, prior=prior2,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
summary(m2)
summary(m3)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(-1, 0, 2,1,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res21-m1 different from expected"))
}
if(any(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m2$Sol, m2$VCV)))[,2]<tpar)/length(tpar), "res21-m2 different from expected"))
}
tpar<-c(-1, 0, 2,0,0,1,1)
if(any(HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m3$Sol, m3$VCV)))[,2]<tpar)/length(tpar), "res21-m3 different from expected"))
}
res21[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
res21b[i,]<-posterior.mode(mcmc(cbind(m2$Sol, m2$VCV)))
res21c[i,]<-posterior.mode(mcmc(cbind(m3$Sol, m3$VCV)))
print(i)
}
psets<-c(psets, c(-1, 0, 2,1,1), c(-1, 0, 2,1,1), c(-1, 0, 2,0,0,1,1))

# caetgorical k=3 + random 5.9 seconds
print("res22")
res22<-matrix(NA, nsim,10)
R=diag(2)
G=matrix(c(2,0.5,0.5,1),2,2)
prior=list(R=list(V=R, n=1, fix=1),G=list(G1=list(V=G, n=1)))
tpar<-c(-1,0,2,0.5,0.5,1,1,0,0,1)

for(i in 1:nsim){
fac<-as.factor(sample(1:100,900,replace=TRUE))
l<-mvrnorm(900, c(-1,0), R)+mvrnorm(100, c(0,0), G)[fac,]

y<-1:900
for(j in 1:900){
y[j]<-which(rmultinom(1,1,prob=c(1,exp(l[j,][1]), exp(l[j,][2])))==1)
}

data=data.frame(y=y, fac=fac)
m1<-MCMCglmm(y~trait-1, random=~us(trait):fac, family=c("categorical"), rcov=~us(trait):units, data=data, prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res22 different from expected"))
}
res22[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

# binary probit

print("res23")
res23<-matrix(NA, nsim,3)
R<-as.matrix(1)
G<-as.matrix(2)
prior=list(R=list(V=R, n=1, fix=1))
tpar<-c(-1,0.5,1)

for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-rnorm(300, -1+0.5*x, R)
cp<-0.5
pr<-cbind(pnorm(-g),1-pnorm(-g))

y1<-1:300
for(j in 1:300){
  y1[j]<-which(rmultinom(1, 1, pr[j,])==1)
}

table(y1)
cp<-qnorm(cumsum(c(0,table(y1)/length(y1))), 0, sqrt(1+R))
cp-cp[2]

data=data.frame(y1=y1, fac=fac, xcov=x)
if(DICtest){
m1<-MCMCglmm(y1~xcov, data=data, prior=prior, family=c("ordinal"), verbose=verbose, nitt=3, thin=1, burnin=1, pl=TRUE, slice=TRUE)
d<-sum(dbinom(data$y1-1, 1, pnorm(m1$Liab[1,]), log=TRUE))
if(abs(-2*d-m1$Deviance[2])<1e-6){
 print("Deviance OK for sliced probit (res23)")
}else{
 stop("Deviance OK for sliced probit (res23)")
}}
m1<-MCMCglmm(y1~xcov,data=data, prior=prior, family="ordinal", verbose=verbose, nitt=nitt, thin=thin, burnin=burnin, slice=TRUE)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res23 different from expected"))
}
res23[i,]<-c(posterior.mode(mcmc(cbind(m1$Sol, m1$VCV))))
print(i)
}
psets<-c(psets, tpar)


# 4 category ordinal

print("res24")
res24<-matrix(NA, nsim,5)
R<-as.matrix(1)
G<-as.matrix(2)
prior=list(R=list(V=R, n=1, fix=1)) 
cp<-c(0.5,1)
tpar<-c(cp,-1,1,1)
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-rnorm(300, -1+1*x, R)
pr<-cbind(pnorm(-g),pnorm(cp[1]-g)-pnorm(-g),  pnorm(cp[2]-g)-pnorm(cp[1]-g), 1-pnorm(cp[2]-g))

y1<-1:300
for(j in 1:300){
  y1[j]<-which(rmultinom(1, 1, pr[j,])==1)
}

data=data.frame(y1=y1, fac=fac, xcov=x)
m1<-MCMCglmm(y1~xcov,data=data, prior=prior, family="ordinal", verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res24 different from expected"))
}
res24[i,]<-c(posterior.mode(mcmc(cbind(m1$CP, m1$Sol, m1$VCV))))
print(i)
}
psets<-c(psets, tpar)

# bivariate ordinal

print("res25")
res25<-matrix(NA, nsim,13)
R<-diag(2)
G<-matrix(c(1,0.5,0.5,1),2,2)
prior=list(R=list(V=R, n=1, fix=1),G=list(G1=list(V=G, n=2)))
for(i in 1:nsim){
fac<-as.factor(sample(1:75,300,replace=TRUE))
x<-runif(300)
g<-mvrnorm(300, c(0,0), R)+mvrnorm(75, c(0,0), G)[fac,]
g[,1]<-g[,1]-1+1*x
g[,2]<-g[,2]-0.5-1*x

cp1<-c(0.5,1)
cp2<-c(0.75)
pr1<-cbind(pnorm(-g[,1]),pnorm(cp1[1]-g[,1])-pnorm(-g[,1]),  pnorm(cp1[2]-g[,1])-pnorm(cp1[1]-g[,1]), 1-pnorm(cp1[2]-g[,1]))
pr2<-cbind(pnorm(-g[,2]),pnorm(cp2[1]-g[,2])-pnorm(-g[,2]),  1-pnorm(cp2[1]-g[,2]))


y1<-1:300
y2<-1:300
for(j in 1:300){
  y1[j]<-which(rmultinom(1, 1, pr1[j,])==1)
  y2[j]<-which(rmultinom(1, 1, pr2[j,])==1)
}

data=data.frame(y1=y1, y2=y2, fac=fac, xcov=x)
m1<-MCMCglmm(cbind(y1,y2)~trait+trait:xcov-1,random=~us(trait):fac, rcov=~idh(trait):units, data=data, prior=prior, family=cbind("ordinal", "ordinal"),verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)), ask=FALSE)
}
tpar<-c(cp1, cp2,-1,-0.5,1,-1,1,0.5,0.5,1,1,1)
if(any(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res25 different from expected"))
}
res25[i,]<-posterior.mode(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

# idv 

print("res26")
res26<-matrix(NA, nsim,3)
R<-diag(1)
G<-diag(1)
prior=list(R=list(V=R, nu=1),G=list(G1=list(V=G, nu=1)))
for(i in 1:nsim){
fac1<-as.factor(sample(1:75,300,replace=TRUE))
fac2<-as.factor(sample(1:75,300,replace=TRUE))
y<-mvrnorm(300, 0, R)+mvrnorm(75, 0, G)[fac1]+mvrnorm(75, 0, G)[fac2]
data=data.frame(y=y, fac1=fac1,fac2=fac2)
m1<-MCMCglmm(y~1,random=~idv(fac1+fac2), data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}

res26[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
tpar<-c(0,1,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res26 different from expected"))
}
psets<-c(psets, tpar)

# idv +iid random effects

print("res27")
res27<-matrix(NA, nsim,6)
R<-diag(1)
G<-diag(1)
G2<-diag(1)
prior=list(R=list(V=R, nu=1),G=list(G1=list(V=G, nu=1, alpha.mu=1, alpha.V=100), G2=list(V=G2, nu=1, alpha.mu=1, alpha.V=100)))

for(i in 1:nsim){
fac1<-as.factor(sample(1:75,300,replace=TRUE))
fac2<-as.factor(sample(1:75,300,replace=TRUE))
id<-as.factor(sample(1:75,300,replace=TRUE))
fac3<-as.factor(sample(1:3, 300, T))
y<-mvrnorm(300, 0, R)+mvrnorm(75, 0, G)[fac1]+mvrnorm(75, 0, G)[fac2]+rnorm(75,0,sqrt(G2))[id]
data=data.frame(y=y, fac1=fac1,fac2=fac2, fac3=fac3, id=id)
m1<-MCMCglmm(y~fac3,random=~idv(fac1+fac2)+id, data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(mcmc(m1$Sol[,1:3]), m1$VCV)), ask=FALSE)
}
tpar<-c(0,0,0,1,1,1)
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res27 different from expected"))
}

res27[i,]<-posterior.mode(mcmc(cbind(m1$Sol[,1:3], m1$VCV)))
print(i)
}
psets<-c(psets, tpar)


# multivariate binary
print("res28")
nT<-2
res28<-matrix(NA, nsim,nT*(1+nT))
R<-cbind(c(1,-0.25),c(-0.25,1))
mu<-matrix(c(0.3,0.9),nT,1)
prior=list(R=list(V=diag(nT), nu=nT+1))

for(i in 1:nsim){
	y<-mvrnorm(300, mu, R)
	data=data.frame(apply(y,2, function(x){rbinom(length(x), 1, plogis(x))}))
	m1<-MCMCglmm(cbind(X1, X2)~trait-1,rcov=~corg(trait):units, family=rep("categorical", nT), data=data, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}	
	if(plotit){
		plot(mcmc(cbind(mcmc(m1$Sol), m1$VCV)), ask=FALSE)
	}
	tpar<-c(mu,c(R))
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res28 different from expected"))
}

	res28[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
	print(i)
}
psets<-c(psets, tpar)

# bivaraite 4 category ordinal
print("res29")
res29<-matrix(NA, nsim,12)
R<-cbind(c(1,-0.35), c(-0.35,1))

prior=list(R=list(V=diag(2), n=3)) #,G=list(G1=list(V=G, n=1)))
for(i in 1:nsim){
	fac<-as.factor(sample(1:75,300,replace=TRUE))
	x<-runif(300)
	g<-mvrnorm(300, c(-1,0), R)+cbind(x*-1, x*0.5) #+mvrnorm(75, c(0), G)[fac]
	cp<-c(0.5,1)
	pr<-cbind(pnorm(-g[,1]),pnorm(cp[1]-g[,1])-pnorm(-g[,1]),  pnorm(cp[2]-g[,1])-pnorm(cp[1]-g[,1]), 1-pnorm(cp[2]-g[,1]))
	cp2<-c(1,2)
	pr2<-cbind(pnorm(-g[,2]),pnorm(cp2[1]-g[,2])-pnorm(-g[,2]),  pnorm(cp2[2]-g[,2])-pnorm(cp2[1]-g[,2]), 1-pnorm(cp2[2]-g[,2]))
	
	y1<-1:300
	for(j in 1:300){
		y1[j]<-which(rmultinom(1, 1, pr[j,])==1)
	}
	y2<-1:300
	for(j in 1:300){
		y2[j]<-which(rmultinom(1, 1, pr2[j,])==1)
	}
	data=data.frame(y1=y1, y2=y2, fac=fac, xcov=x)
	m1<-MCMCglmm(cbind(y1, y2)~trait+trait:xcov-1,rcov=~corg(trait):units, data=data, prior=prior, family=c("ordinal", "ordinal"), verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
	if(plotit){
		plot(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)), ask=FALSE)
	}
	tpar<-c(cp,cp2, -1,0,-1,0.5, c(R))
        if(any(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$CP, m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res29 different from expected"))
        }
	res29[i,]<-c(posterior.mode(m1$CP),posterior.mode(m1$Sol),posterior.mode(m1$VCV))
	print(i)
}
psets<-c(psets, tpar)

# hurdle Poisson

print("res30")
res30<-matrix(NA, nsim,6)

tpar<-c(-1,log(1/2),1,1/2,1,1)

for(i in 1:nsim){

  x<-rnorm(300)
  l1<-rnorm(300, -1+x, sqrt(1))
  l2<-rnorm(300, log(1/2)+x/2, sqrt(1))
     
  y<-rbinom(300, 1, 1-plogis(l2))
  y[which(y==1)]<-qpois(runif(sum(y==1), dpois(0, exp(l1[which(y==1)])), 1), exp(l1[which(y==1)]))  
  # cunning sampler from Peter Dalgaard (R-sig-mixed)


data=data.frame(y=y, x=x)
prior=list(R=list(V=diag(2), fix=2, nu=1))
m1<-MCMCglmm(y~trait-1+trait:x, rcov=~idh(trait):units, data=data, family="hupoisson", prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

res30[i,]<-c(posterior.mode(m1$Sol), posterior.mode(m1$VCV))
if(SUMtest){
summary(m1)
}
	if(plotit){
		plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
	}
        if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res30 different from expected"))
        }
print(i)
}
psets<-c(psets, tpar)


# truncated Poisson

print("res31")
res31<-matrix(NA, nsim,3)

tpar<-c(-1,1/2, 1)
for(i in 1:nsim){

 
 l<-rnorm(300, -1+x/2, 1)
 y<-qpois(runif(300, dpois(0, exp(l)), 1), exp(l))  

dat<-data.frame(y=y,x=x)
m1<-MCMCglmm(y~x, family="ztpoisson", data=dat, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

res31[i,]<-c(posterior.mode(m1$Sol),posterior.mode(m1$VCV))

        if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,2]<tpar)/length(tpar), "res31 different from expected"))
        }
print(i)
}
psets<-c(psets, tpar)


# geometric

print("res32")
res32<-matrix(NA, nsim,2)
tpar<-c(-1,0.5)
for(i in 1:nsim){


nu<--1
v<-0.5

y<-rgeom(300,plogis(rnorm(300,nu,sqrt(v))))


dat<-data.frame(y=y)
m1<-MCMCglmm(y~1, family="geometric", data=dat, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
res32[i,]<-c(posterior.mode(m1$Sol),posterior.mode(m1$VCV))

        if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV[,1])))[,2]<tpar)/length(tpar), "res32 different from expected"))
        }
print(i)
}
psets<-c(psets, tpar)

#zibinomial

print("res33")
res33<-matrix(0, nsim, 6)
R<-diag(2)
prior=list(R=list(V=R, n=1, fix=2))
tune=list(diag(2))
tpar<-c(0, -1, 1, 1/2, 1, 1)
for(i in 1:nsim){
	
x<-rnorm(300)	
l<-mvrnorm(300,c(0,-1), R)
l[,1]<-l[,1]+x
l[,2]<-l[,2]+x/2
y<-rbinom(300, 1, 1-plogis(l[,2]))
y[which(y==1)]<-rbinom(sum(y==1), 20, plogis(l[,1][which(y==1)]))  


data=data.frame(success=y,failure=20-y, x=x)
m1<-MCMCglmm(cbind(success,failure)~trait-1+trait:x, rcov=~idh(trait):units, data=data, family="zibinomial",prior=prior,verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res33 different from expected"))
}
res33[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

print("res34")
# sir model
res34<-matrix(0, nsim, 4)
tpar<-c(0,0,-0.5, 1)
for(i in 1:nsim){

id <- sample(1:100, 100, replace = T)
y <- rnorm(100)

id1 <- sample(id, 100, replace = T)
id <- as.factor(id)
id1 <- factor(id1, levels=levels(id))

L<-diag(100)-sir(~id1, ~id)*-0.5

y<-solve(L,y)

my.data <- data.frame(y = y, id = id, id1 = id1, x = rnorm(100))

if(DICtest){
  m1 <- MCMCglmm(y ~ x + sir(~id1, ~id), data = my.data, verbose=verbose, nitt=3, thin=1, burnin=0)

  L<-Diagonal(100)-m1$XL*m1$Lambda[2]

  d<-sum(dnorm((L%*%y)@x, (m1$X%*%m1$Sol[2,])@x, sqrt(m1$VCV[2]), log=TRUE))+log(det(L))

  if(abs(-2*d-m1$Deviance[2])<1e-6){
    print("Deviance OK for sir model (res34)")
  }else{
    stop("Deviance wrong  sir model (res34)")
  }
}


m1 <- MCMCglmm(y ~ x + sir(~id1, ~id), data = my.data, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
res34[i,]<-c(posterior.mode(m1$Sol),posterior.mode(m1$Lambda), posterior.mode(m1$VCV))
        if(any(HPDinterval(mcmc(cbind(m1$Sol,m1$Lambda, m1$VCV[,1])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$Lambda,m1$VCV[,1])))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$Lambda, m1$VCV[,1])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol,m1$Lambda, m1$VCV[,1])))[,2]<tpar)/length(tpar), "res34 different from expected"))
        }


print(i)
}

psets<-c(psets, tpar)


 # Univariate threshold probit
print("res35") 

res35<-matrix(0, nsim, 3)

prior<-list(R=list(V=as.matrix(1), n=1, fix=1))
tpar<-c(0.5,-1, 1)
for(i in 1:nsim){

prior<-list(R=list(V=as.matrix(1), n=1, fix=1))

x<-rnorm(100)
y<-rnorm(100, 0.5-x)
data<-data.frame(y=as.numeric(y>0), x=x)
cp<-c(-Inf, 0, Inf)
if(DICtest){
m1<-MCMCglmm(y~x, family="threshold",  data=data, prior=prior, verbose=verbose, nitt=3, thin=1, burnin=1)
d<-sum(log(pnorm(cp[data$y+2], (m1$X%*%m1$Sol[2,])@x, sqrt(m1$VCV[2,1]))-pnorm(cp[data$y+1], (m1$X%*%m1$Sol[2,])@x, sqrt(m1$VCV[2,1]))))
if(abs(-2*d-m1$Deviance[2])<1e-6){
 print("Deviance OK for univariate threshold (res35)")
}else{
 stop("Deviance wrong for univariate threshold (res35)")
}}
m1<-MCMCglmm(y~x, family="threshold", data=data, verbose=verbose, prior=prior, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res35 different from expected"))
}
res35[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)


print("res36")  # Test threshold against probit

res36<-matrix(0, nsim, 8)

R<-cbind(c(2,1), c(1,1))
tpar<-c(1,0.5, -1, -0.5, c(R))

prior<-list(R=list(V=R, nu=0.002, fix=2))

for(i in 1:nsim){

y<-mvrnorm(100, c(1, 0.5), R)
x<-rnorm(100)
y[,1]<-y[,1]-x
y[,2]<-y[,2]-0.5*x

data=data.frame(y1=y[,1], y2=as.numeric(y[,2]>0), x=x)

m1<-MCMCglmm(cbind(y1,y2)~trait-1+trait:x, rcov=~us(trait):units, family=c("gaussian", "threshold"), data=data, verbose=verbose, prior=prior, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res36 different from expected"))
}
res36[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

print("res37")  # Test threshold against probit

res37<-matrix(0, nsim, 3)

prior<-list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))
tpar<-c(0,1,1)

for(i in 1:nsim){
tree<-rcoal(100)
y<-rbv(tree, diag(1), nodes="TIPS")+rnorm(100)
A<-inverseA(tree)

data=data.frame(y=as.numeric(y>0), species=tree$tip.label)

m1<-MCMCglmm(y~1, random=~species, rcov=~units, family="threshold", ginv=list(species=A$Ainv),verbose=verbose,  data=data, prior=prior, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res37 different from expected"))
}
res37[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

res38<-matrix(0, nsim, 12)

tpar<-c(1,0,0.5, -0.5, 1,-0.5,-0.5,1, 1,0.25,0.25,1)

G<-matrix(c(1,-0.5,-0.5,1),2,2)
R<-matrix(c(1, 0.25, 0.25,1),2,2)

for(i in 1:nsim){
id<-c(1:50, rep(51:100,4))
x<-rnorm(250)
y<-mvrnorm(100, c(0,0), G)[id,]+mvrnorm(250, c(0,0), R)
y[,1]<-y[,1]+1+0.5*x
y[,2]<-y[,2]+0-0.5*x

data=data.frame(y1=as.numeric(y[,1]>0), y2=as.numeric(y[,2]>0), id=id, x=x)

prior=list(R=list(V=diag(2), nu=10), G=list(G1=list(V=G, nu=3)))

m1<-MCMCglmm(cbind(y1, y2)~trait-1+trait:x, random=~us(trait):id, rcov=~corg(trait):units, family=c("threshold","threshold"), data=data, verbose=verbose, prior=prior, nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res38 different from expected"))
}
res38[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

res39<-matrix(0, nsim, 12)
print("res39")  # Block-diagonal R-structure with augmented animal model
data(BTped)

tpar<-c(0,0,1, 0, 0,1,1,0, 0,1,1,2)

prior<-list(G=list(G1=list(V=diag(2), nu=1.002)), R=list(R1=list(V=diag(2), nu=1.002), R2=list(V=diag(2), nu=1.002)))

for(i in 1:nsim){

fac<-as.factor(sample(1:2, sum(!is.na(BTped[,2])), TRUE))

y1<-rnorm(sum(!is.na(BTped[,2])), 0, sqrt(as.numeric(fac)))
y1<-y1+rbv(BTped, 1)[which(!is.na(BTped[,2]))]

y2<-rnorm(sum(!is.na(BTped[,2])), 0, sqrt(as.numeric(fac)))
y2<-y2+rbv(BTped, 1)[which(!is.na(BTped[,2]))]

dat<-data.frame(y1=y1, y2=y2, fac=fac, id=BTped[,1][which(!is.na(BTped[,2]))])

A<-inverseA(BTped)

m1<-MCMCglmm(cbind(y1,y2)~trait-1, random=~us(trait):id, rcov=~us(trait:at.level(fac,1)):units+idh(trait:at.level(fac,2)):units, data=dat, family=rep("gaussian", 2), prior=prior, verbose=verbose,nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res39 different from expected"))
}
res39[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

res40<-matrix(0, nsim, 2)
print("res40")  # Reduced phylogenetic probit model
data(BTped)

tpar<-c(0.5, -1)

   data(bird.families) 
   A<-inverseA(bird.families, reduced=TRUE)
   dat<-data.frame(taxon=bird.families$tip.label)
   dat$mii<-A$mii[match(dat$taxon,A$pedigree[,1])]
   dat$mother<-A$pedigree[,2][match(dat$taxon,A$pedigree[,1])]

for(i in 1:nsim){

    
     l<-0.5+rbv(bird.families, 1, nodes="TIPS") 
     x<-rnorm(length(l))     
     y<-(l-x)>0
          
     dat$y<-y
     dat$x<-x

     prior<-list(R=list(V=diag(dat$mii), fix=1), G=list(G1=list(V=1, fix=1)))

     # inverse matrix of shared phyloegnetic history
          
     m1<-MCMCglmm(y~x, random=~mother, rcov=~idh(units):units, ginverse=list(mother=A$Ainv),
     data=dat, prior=prior,  verbose=verbose, nitt=nitt, thin=thin, burnin=burnin, family="threshold")



if(SUMtest){
summary(m1)
}
if(plotit){
plot(mcmc(m1$Sol), ask=FALSE)
}
if(any(HPDinterval(mcmc(m1$Sol))[,1]>tpar | HPDinterval(mcmc(m1$Sol))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(m1$Sol))[,1]>tpar | HPDinterval(mcmc(m1$Sol))[,2]<tpar)/length(tpar), "res40 different from expected"))
}
res40[i,]<-posterior.mode(mcmc(m1$Sol))
print(i)
}
psets<-c(psets, tpar)


print("res41")  # Trivariate Gaussian, Poisson and 3 bin threshold

res41<-matrix(0, nsim, 16)

R<-diag(c(2,1,1))
R[1,3]<-R[3,1]<-1
cp<-1
tpar<-c(1,0,0.5, -1, 0,-0.5, c(R), cp)

prior<-list(R=list(V=R, nu=0.002, fix=3))

for(i in 1:nsim){

y<-mvrnorm(100, c(1, 0, 0.5), R)
x<-rnorm(100)
y[,1]<-y[,1]-x
y[,2]<-y[,2]
y[,3]<-y[,3]-0.5*x


data=data.frame(y1=y[,1], y2=rpois(100, exp(y[,2])), y3=as.numeric(cut(y[,3], c(-Inf, 0, cp, Inf))), x=x)
if(DICtest){
m1<-MCMCglmm(cbind(y1,y2, y3)~trait-1+trait:x, rcov=~us(trait):units, family=c("gaussian", "poisson", "threshold"), data=data, verbose=verbose, prior=prior, nitt=3, thin=1, burnin=1,  pl=TRUE)

pred1<-(m1$X%*%m1$Sol[2,])[1:100]
pred2<-(m1$X%*%m1$Sol[2,])[101:200]
pred3<-(m1$X%*%m1$Sol[2,])[201:300]

l1<-m1$Liab[1,1:100]
l2<-m1$Liab[1,101:200]
l3<-m1$Liab[1,201:300]

tcp<-c(-Inf, 0, m1$CP[2], Inf)

Ritt<-matrix(m1$VCV[2,],3,3)
cpred3<-pred3+t((Ritt[3,1:2]%*%solve(Ritt[1:2,1:2]))%*%t(cbind(c(l1-pred1), c(l2-pred2))))
cR3<-Ritt[3,3]-Ritt[3,1:2]%*%solve(Ritt[1:2,1:2])%*%Ritt[1:2,3]

cpred1<-pred1+Ritt[1,2]%*%solve(Ritt[2,2])%*%(l2-pred2)
cR1<-Ritt[1,1]-Ritt[1,2]%*%solve(Ritt[2,2])%*%Ritt[1,2]

d<-sum(log(pnorm(tcp[data$y3+1], cpred3, sqrt(cR3))-pnorm(tcp[data$y3], cpred3, sqrt(cR3))))
d<-d+sum(dpois(data$y2, exp(l2), log=TRUE))
d<-d+sum(dnorm(data$y1, cpred1, sqrt(cR1), log=TRUE))

if(abs(-2*d-m1$Deviance[2])<1e-6){
 print("Deviance OK for trivariate threshold/gaussian/non-gaussian (res41)")
}else{
 stop("Deviance OK for trivariate threshold/gaussian/non-gaussian (res41)")
}}

m1<-MCMCglmm(cbind(y1,y2, y3)~trait-1+trait:x, rcov=~us(trait):units, family=c("gaussian", "poisson", "threshold"), data=data, verbose=verbose, prior=prior, nitt=nitt, thin=thin, burnin=burnin)
if(SUMtest){
summary(m1)
}

if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV, m1$CP)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV, m1$CP)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV, m1$CP)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV, m1$CP)))[,2]<tpar)/length(tpar), "res41 different from expected"))
}
res41[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV, m1$CP)))
print(i)
}
psets<-c(psets, tpar)

print("res42")  # Trivariate: 1 Gaussian, 2 threshold with correlation sub-matrix

res42<-matrix(0, nsim, 15)

R<-diag(c(2,1,1))
R[1,3]<-R[3,1]<-1
R[2,3]<-R[3,2]<-0.25

tpar<-c(1,0,0.5, -1, 0,-0.5, c(R))

prior<-list(R=list(V=R, nu=0, fix=2))

for(i in 1:nsim){

y<-mvrnorm(100, c(1, 0, 0.5), R)
x<-rnorm(100)
y[,1]<-y[,1]-x
y[,2]<-y[,2]
y[,3]<-y[,3]-0.5*x


data=data.frame(y1=y[,1], y2=y[,2]>0, y3=y[,3]>0, x=x)

if(DICtest){print("DIC not available for this model")}
m1<-MCMCglmm(cbind(y1,y2, y3)~trait-1+trait:x, rcov=~cors(trait):units, family=c("gaussian", "threshold", "threshold"), data=data, prior=prior, verbose=verbose,nitt=nitt, thin=thin, burnin=burnin)


if(SUMtest){
summary(m1)
}

if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res42 different from expected"))
}
res42[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

print("res43")  # 1st order ante-dependence on 4 observations

res43<-matrix(0, nsim, 9)

R<-diag(4:1)
beta<-c(-0.5,0,0.5)

tpar<-c(1,0.5, diag(R), beta)

Lambda1<-Diagonal(4,1)-bandSparse(4,4,-1,list(beta))
V<-solve(Lambda1)%*%R%*%t(solve(Lambda1))


for(i in 1:nsim){

y<-mvrnorm(100, rep(1,4), V)
x<-rnorm(100*4)
y<-y+x*0.5

data=data.frame(y=c(y), day=as.factor(rep(1:4,each=100)), id=as.factor(rep(1:100,4)), x=x)


if(DICtest){
m1<-MCMCglmm(y~x, rcov=~ante1(day):id, data=data, verbose=verbose, nitt=3, thin=1, burnin=1)
dev<-0
for(j in 1:100){
dev<-dev+mvtnorm::dmvnorm(data$y[which(data$id==j)], m1$Sol[1,1]+m1$Sol[1,2]*data$x[which(data$id==j)], matrix(m1$VCV[1,],4,4), log=TRUE)
}
if(abs(-2*dev-m1$Deviance[1])<1e-6){
 print("Deviance OK for ante-dependence Gaussian (res43)")
}else{
 stop("Deviance wrong for ante-dependence Gaussian (res43)")
}}

m1<-MCMCglmm(y~x, rcov=~ante1(day):id, data=data, verbose=verbose,nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}

if(plotit){
plot(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV))), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV))))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV))))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV))))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV))))[,2]<tpar)/length(tpar), "res43 different from expected"))
}
res43[i,]<-posterior.mode(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV))))
print(i)
}
psets<-c(psets, tpar)

print("res44")  # 2nd order ante-dependence with constant innovation variance on 4 observations

res44<-matrix(0, nsim, 8)

R<-diag(4)
beta<-c(-0.5,0,0.5,0,0)

tpar<-c(1,0.5, 1, beta)

Lambda1<-Diagonal(4,1)-bandSparse(4,4,-(1:2),list(beta[1:3], beta[4:5]))
V<-solve(Lambda1)%*%R%*%t(solve(Lambda1))


for(i in 1:nsim){

y<-mvrnorm(100, rep(1,4), V)
x<-rnorm(100*4)
y<-y+x*0.5

data=data.frame(y=c(y), day=as.factor(rep(1:4,each=100)), id=as.factor(rep(1:100,4)), x=x)

m1<-MCMCglmm(y~x, rcov=~ante2v(day):id, data=data, verbose=verbose,nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}

if(plotit){
plot(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,2)[,-c(2:4)])), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,2)[,-c(2:4)])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,2)[,-c(2:4)])))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,2)[,-c(2:4)])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,2)[,-c(2:4)])))[,2]<tpar)/length(tpar), "res44 different from expected"))
}
res44[i,]<-posterior.mode(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,2)[,-c(2:4)])))
print(i)
}
psets<-c(psets, tpar)

print("res45")  # 3rd order ante-dependence with constant innovation variance and regressions on 4 observations

res45<-matrix(0, nsim, 6)

R<-diag(4)
beta<-c(-0.5,0,0.5)

tpar<-c(1,0.5, 1, beta)

Lambda1<-Diagonal(4,1)-bandSparse(4,4,-(1:3),list(rep(beta[1],3), rep(beta[2],2), beta[3]))

V<-solve(Lambda1)%*%R%*%t(solve(Lambda1))


for(i in 1:nsim){

y<-mvrnorm(100, rep(1,4), V)
x<-rnorm(100*4)
y<-y+x*0.5

data=data.frame(y=c(y), day=as.factor(rep(1:4,each=100)), id=as.factor(rep(1:100,4)), x=x)

m1<-MCMCglmm(y~x, rcov=~antec3v(day):id, data=data, verbose=verbose,nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}

if(plotit){
plot(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,3)[,c(1,5,8,10)])), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,3)[,c(1,5,8,10)])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,3)[,c(1,5,8,10)])))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,3)[,c(1,5,8,10)])))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,3)[,c(1,5,8,10)])))[,2]<tpar)/length(tpar), "res44 different from expected"))
}
res45[i,]<-posterior.mode(mcmc(cbind(m1$Sol, posterior.ante(m1$VCV,3)[,c(1,5,8,10)])))
print(i)
}
psets<-c(psets, tpar)


print("res46")  # 3-category multinomial

res46<-matrix(0, nsim, 6)

I<-diag(2)
J<-matrix(1,2,2)
R<-(I+J)/3

tpar<-c(-1,-2,c(R))

for(i in 1:nsim){

y<-mvrnorm(300, c(-1,1), R)

prob<-cbind(1,exp(y))/(1+rowSums(exp(y)))

y<-t(apply(prob, 1, function(x){rmultinom(1, size=sample(1:10), prob=x)}))

data=data.frame(y1=y[,1], y2=y[,2], y3=y[,3])

m1<-MCMCglmm(cbind(y1,y2,y3)~trait-1, rcov=~us(trait):units, data=data, family="multinomial3", verbose=verbose,nitt=nitt, thin=thin, burnin=burnin)

if(SUMtest){
summary(m1)
}

if(plotit){
plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
}
if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res46 different from expected"))
}
res46[i,]<-posterior.mode(mcmc(cbind(m1$Sol, m1$VCV)))
print(i)
}
psets<-c(psets, tpar)

# zero-altered Poisson

print("res47")
res47<-matrix(NA, nsim,5)

tpar<-c(1,0,1, 0, 1)

for(i in 1:nsim){

  x<-rnorm(300)
  l1<-rnorm(300, 1+x, sqrt(1))
  l2<-rnorm(300, 1+x, sqrt(1))

  y<-rbinom(300, 1, 1-exp(-exp(l2)))
  y[which(y==1)]<-qpois(runif(sum(y==1), dpois(0, exp(l1[which(y==1)])), 1), exp(l1[which(y==1)]))  
  # cunning sampler from Peter Dalgaard (R-sig-mixed)

data=data.frame(y=y, x=x)
prior=list(R=list(V=diag(1), nu=1))
m1<-MCMCglmm(y~trait*x, rcov=~trait:units, data=data, family="zapoisson", prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)

res47[i,]<-c(posterior.mode(m1$Sol), posterior.mode(m1$VCV))
if(SUMtest){
summary(m1)
}
	if(plotit){
		plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
	}
        if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res47 different from expected"))
        }
print(i)
}
psets<-c(psets, tpar)


# exponential

print("res48")
res48<-matrix(NA, nsim,2)

tpar<-c(-1,1)

for(i in 1:nsim){



  l<-rnorm(300, -1, sqrt(1))
     
  y<-rexp(300, rate=exp(l))



data=data.frame(y=y)
prior=list(R=list(V=diag(1), nu=1))
m1<-MCMCglmm(y~1,data=data, family="exponential", prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
res48[i,]<-c(posterior.mode(m1$Sol), posterior.mode(m1$VCV))
if(SUMtest){
summary(m1)
}
	if(plotit){
		plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
	}
        if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res48 different from expected"))
        }
print(i)
}
psets<-c(psets, tpar)


# covu

print("res49")
res49<-matrix(NA, nsim,13)

tpar<-c(1,1,-1, 1,0.5, 0.25, 0.5,1, -0.25, 0.25, -0.25, 0.5,1)

for(i in 1:nsim){

  n<-100
  V<-matrix(c(1,0.5, 0.25, 0.5,1, -0.25, 0.25, -0.25, 0.5),3,3)
  Vr<-1

  u<-mvrnorm(n, c(0,0,0), V)

  ya<-1+u[,2]
  yb<-1+u[,3]
  individual<-as.factor(rep(1:n, 4))
  type<-as.factor(c(rep("s", 2*n), rep("r",2*n)))
  measure<-as.factor(c(rep("a", n),rep("b", n), rep("c",2*n)))
  yc<--1+u[individual[which(measure=="c")],1]+rnorm(2*n,0,sqrt(Vr))

  dat<-data.frame(y=c(ya,yb,yc), type=type, individual=individual, measure=measure)
 
  prior<-list(R=list(R1=list(V=V, nu=3, covu=TRUE), R2=list(V=Vr, nu=1)))

    if(DICtest){
    m1<-MCMCglmm(y~measure-1, random=~us(at.level(type,"r")):individual, rcov=~us(at.level(type, "s"):measure):individual+us(at.level(type, "r")):units, data=dat, prior=prior, pr=TRUE, nitt=2, thin=1, burnin=1, verbose=FALSE)

   Vest<-matrix(m1$VCV[1,1:9],3,3)

   Vreg<-Vest[2:3,1]%*%solve(Vest[1,1])
   Vres<-Vest[2:3,2:3]-Vest[2:3,1]%*%solve(Vest[1,1])%*%Vest[1,2:3]
   dev<-0

   for(j in 1:n){
     dev<-dev+mvtnorm::dmvnorm(cbind(ya,yb)[j,], m1$Sol[1,1:2]+Vreg%*%m1$Sol[1,3+j], Vres, log=TRUE)
   }
   dev<-dev+sum(dnorm(yc, m1$Sol[1,3]+m1$Sol[1,3+1:n][individual[which(measure=="c")]], sqrt(m1$VCV[1,10]), log=TRUE))

   if(abs(-2*dev-m1$Deviance[1])<1e-6){
     print("Deviance OK for covu (res49)")
   }else{
     stop("Deviance wrong for covu (res49)")
   }}

   m1<-MCMCglmm(y~measure-1, random=~us(at.level(type,"r")):individual, rcov=~us(at.level(type, "s"):measure):individual+us(at.level(type, "r")):units, data=dat, prior=prior, verbose=verbose, nitt=nitt, thin=thin, burnin=burnin)
   if(SUMtest){
    summary(m1)
   }
   res49[i,]<-c(posterior.mode(m1$Sol), posterior.mode(m1$VCV))

	if(plotit){
		plot(mcmc(cbind(m1$Sol, m1$VCV)), ask=FALSE)
	}
        if(any(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)){
        print(paste(sum(HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,1]>tpar | HPDinterval(mcmc(cbind(m1$Sol, m1$VCV)))[,2]<tpar)/length(tpar), "res49 different from expected"))
        }
print(i)

 }
psets<-c(psets, tpar)






est<-colMeans(cbind(res1, res2, res3, res3b, res4, res4c, res5,res5b,res6, res7,res7b, res8,res9, res10, res11, res12, res13, res14, res15, res17a,res17b, res18, res19, res19b, res20, res21, res21b, res21c, res22, res23, res24, res25, res26, res27, res28, res29, res30, res31, res32, res33, res34, res35, res36, res37, res38, res39, res40, res41, res42, res43, res44, res45, res46, res47, res48, res49), na.rm=T)

np<-c(ncol(res1), ncol(res2), ncol(res3), ncol(res3b), ncol(res4), ncol(res4c), ncol(res5),ncol(res5b),ncol(res6), ncol(res7),ncol(res7b), ncol(res8),ncol(res9), ncol(res10), ncol(res11), ncol(res12), ncol(res13), ncol(res14), ncol(res15), ncol(res17a),ncol(res17b), ncol(res18), ncol(res19), ncol(res19b), ncol(res20), ncol(res21), ncol(res21b), ncol(res21c), ncol(res22), ncol(res23), ncol(res24), ncol(res25), ncol(res26), ncol(res27), ncol(res28), ncol(res29), ncol(res30), ncol(res31), ncol(res32), ncol(res33), ncol(res34), ncol(res35), ncol(res36), ncol(res37), ncol(res38), ncol(res39), ncol(res40), ncol(res41), ncol(res42), ncol(res43), ncol(res44), ncol(res45), ncol(res46), ncol(res47), ncol(res48), ncol(res49))

nam<-c("res1", "res2", "res3", "res3b", "res4", "res4c", "res5","res5b","res6", "res7","res7b", "res8","res9", "res10", "res11", "res12", "res13", "res14", "res15", "res17a","res17b", "res18", "res19", "res19b", "res20", "res21", "res21b", "res21c", "res22", "res23", "res24", "res25", "res26", "res27", "res28", "res29", "res30", "res31", "res32", "res33", "res34", "res35", "res36", "res37", "res38", "res39", "res40", "res41", "res42", "res43", "res44", "res45", "res46", "res47", "res48", "res49")

nam<-paste(rep(nam,np), unlist(sapply(np,function(x){1:x})), sep=".")

names(est)<-nam

plot(est~psets)
abline(0,1)



