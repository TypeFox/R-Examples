### R code from vignette source 'survival.Rnw'

###################################################
### code chunk number 1: survival.Rnw:76-77
###################################################
options(width=70) 


###################################################
### code chunk number 2: survival.Rnw:79-84
###################################################
library(survival)
T=rexp(25)
delta=sample(0:1,25,TRUE)
y=Surv(T,delta)
y


###################################################
### code chunk number 3: survival.Rnw:123-198
###################################################

library(VGAM)
library(partDSA)
set.seed(1)
p=5
tr.n=250
ts.n=5000
C_level=0.3

x=matrix(NA,tr.n,p)
for (j in 1:p){
x[,j]=sample(1:100,tr.n,replace=TRUE)
}
x=data.frame(x)
names(x)=c(paste("X",1:p,sep=""))

x.test=matrix(NA,ts.n,p)
for (j in 1:p){
x.test[,j]=sample(1:100,ts.n,replace=TRUE)
}
x.test=data.frame(x.test)
names(x.test)=c(paste("X",1:p,sep=""))

scale=ifelse(x[,1]>50|x[,2]>75,5,.5)
time=data.frame(rgpd(tr.n,0,scale,0))
scale.test=ifelse(x.test[,1]>50|x.test[,2]>75,5,.5)
time.test=(rgpd(ts.n,0,scale.test,0))

group2=which(x[,1]>50|x[,2]>75)
group1=c(1:tr.n)[-group2]

level1=0
level2=0
cens.time=NULL
while((abs(C_level-level1)>.02)|(abs(C_level-level2)>.02)){
 cens1.time=runif(length(group1),0,1.25)
   cens2.time=runif(length(group2),0,15)
   
   cens1=apply(cbind(cens1.time,time[group1,]),1,which.min)-1
   cens2=apply(cbind(cens2.time,time[group2,]),1,which.min)-1
  
   level1=1-sum(cens1)/length(cens1)
   level2=1-sum(cens2)/length(cens2)
  
}
cens.time[group1]=cens1.time
cens.time[group2]=cens2.time
y0=apply(cbind(cens.time,time),1,min)
cens0=apply(cbind(cens.time,time),1,which.min)-1

L = quantile(y0,.95)
y = pmin(y0,L)
cens = cens0
cens[y0 > y] = 1

y.new=y
y=log(y)


y.new.test <- time.test
y.new.test <- pmin(y.new.test,L)
y.test=log(y.new.test)
cens.test=rep(1,ts.n)


wt=rep(1,tr.n)
wt.test=rep(1,ts.n)
detach(package:VGAM)

model.KM.IPCW=partDSA(x=x,y=Surv(y,cens),x.test=x.test,y.test=Surv(y.test,cens.test),control=DSA.control(vfold=5, MPD=.05, minsplit=40,minbuck=15, loss.function="IPCW",wt.method="KM",missing="no",cut.off.growth=4))
model.Cox.IPCW=partDSA(x=x,y=Surv(y,cens),x.test=x.test,y.test=Surv(y.test,cens.test),control=DSA.control(vfold=5, MPD=.05, minsplit=40, minbuck=15, loss.function="IPCW",wt.method="Cox",missing="no",cut.off.growth=4))
 





###################################################
### code chunk number 4: survival.Rnw:213-216
###################################################
model.KM.IPCW

model.Cox.IPCW


###################################################
### code chunk number 5: survival.Rnw:243-337
###################################################

tr.n=250
ts.n=5000
p=5
C_level=0.3

set.seed(1)

library(VGAM)

x=matrix(NA,tr.n,p)
for (j in 1:p){
x[,j]=runif(tr.n)
}
x=data.frame(x)
names(x)=c(paste("X",1:p,sep=""))

x.test=matrix(NA,ts.n,p)
for (j in 1:p){
x.test[,j]=runif(ts.n)
}
x.test=data.frame(x.test)
names(x.test)=c(paste("X",1:p,sep=""))

theta=4*as.numeric(x[,1]<=.5 | x[,2]>.5)
psi=exp(theta)
shape=0
scale=1/psi
location=0
time=data.frame(rgpd(tr.n,location, scale,shape))
group2=which(scale==unique(scale)[1])
group1=which(scale==unique(scale)[2])

theta=4*as.numeric(x.test[,1]<=.5 | x.test[,2]>.5)
psi=exp(theta)
shape=0
scale.test=1/psi
location=0
time.test=rgpd(ts.n,location, scale.test,shape)
level1=0
level2=0
cens.time=NULL
upper1=3.5
upper2=.1

while((abs(C_level-level1)>.02)|(abs(C_level-level2)>.02)){
   
   cens1.time=runif(length(group1),0,upper1)
   cens2.time=runif(length(group2),0,upper2)
   
   cens1=apply(cbind(cens1.time,time[group1,]),1,which.min)-1
   cens2=apply(cbind(cens2.time,time[group2,]),1,which.min)-1
  
   level1=1-sum(cens1)/length(cens1)
   level2=1-sum(cens2)/length(cens2)
   if(abs(level1-C_level)>.1){ upper1=ifelse(level1>C_level,upper1+.1,upper1-.05)}
   if(abs(level2-C_level)>.1){upper2=ifelse(level2>C_level,upper2+.01,upper2-.01)}
 
}


if(C_level==0){
y0=time[,1]
cens0=rep(1,tr.n)
}else{
cens.time[group1]=cens1.time
cens.time[group2]=cens2.time
y0=apply(cbind(cens.time,time),1,min)
cens0=apply(cbind(cens.time,time),1,which.min)-1
}

L = quantile(y0,.95)
y = pmin(y0,L)
cens = cens0
cens[y0 > y] = 1


y.new=y
y=log(y)


y.new.test <- time.test
y.test=log(time.test)
cens.test=rep(1,ts.n)


wt=rep(1,tr.n)
wt.test=rep(1,ts.n)
detach(package:VGAM)


model1.Brier=partDSA(x=x,y=Surv(y,cens),x.test=x.test,y.test=Surv(y.test,cens.test),control=DSA.control(vfold=5, MPD=.05, minsplit=40, minbuck=15, loss.function="Brier",missing="no",cut.off.growth=4,brier.vec=-3.9237))
model2.Brier=partDSA(x=x,y=Surv(y,cens),x.test=x.test,y.test=Surv(y.test,cens.test),control=DSA.control(vfold=5, MPD=.05, minsplit=40, minbuck=15, loss.function="Brier",missing="no",cut.off.growth=4,brier.vec=quantile(y,c(.25,.5,.75))))
 


###################################################
### code chunk number 6: survival.Rnw:350-353
###################################################
model1.Brier

model2.Brier


###################################################
### code chunk number 7: survival.Rnw:365-387
###################################################
data("GBSG2", package = "TH.data")
set.seed(1)
data(GBSG2)
dataset=GBSG2
y.new=dataset$time
y.new=ifelse(dataset$time>2000,2000,dataset$time)
y=log(y.new)
cens=dataset$cens
cens[which(dataset$time>2000)]=1
x=dataset[,c(1:8)]
brier.vec=median(y)

set.seed(1)
model.IPCW=partDSA(x=x,y=Surv(y,cens),control=DSA.control(vfold=5, cut.off.growth=5,MPD=.01,minsplit=50,minbuck=20,loss.function="IPCW",wt.method="KM"))

set.seed(1)
model.Brier=partDSA(x=x,y=Surv(y,cens),control=DSA.control(vfold=5, cut.off.growth=5,MPD=.01,minsplit=50,minbuck=20,loss.function="Brier",brier.vec=brier.vec))


model.IPCW

model.Brier


