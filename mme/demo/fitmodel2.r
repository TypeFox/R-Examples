# THIS IS A COMPLETE EXAMPLE WITH MODEL 2
# USE SIMULATED DATA FOR 50 AREAS AND 6 PERIODS OF TIME

######################
# Model 2
# Model with two independent random effects in each category of the response variable:
# one domain random effect and another independent time and domain random effect
######################

library(mme)

#THIS FUNCTION SIMULATE THE DATA
simulation<-function(d,t,k){
D=d*t
u=matrix(0,d,t)
x1=matrix(0,d,t)
x2=matrix(0,d,t)
u1=matrix(0,d,t)
u2=matrix(0,d,t)
for (i in 1:d){
	for (j in 1:t){
		u1[i,j]=((i-d)/d+1/2+j/t)/3
		u2[i,j]=((i-d)/d+2/2+j/t)/3
		x1[i,j]=1+u1[i,j]
		x2[i,j]=1+sqrt(2)*(0*u1[i,j]+sqrt(1-(0*0))*u2[i,j])
}}
phi1=c(1,2)
phi2=c(0.25,0.50)
u1=matrix(0,d,k-1)
s = 12345678
set.seed(s)
u1[,1]=rnorm(d,mean=0,sd=sqrt(phi1[1]))
u1[,2]=rnorm(d,mean=0,sd=sqrt(phi1[2]))

u2=matrix(0,D,k-1)
rho=c(0.50,0.75)
a=omega(t,k,rho,phi2)
ceros=matrix(rep(0,t),t,1)
datos=mvrnorm(d,ceros,((phi2[1])*(a[[1]][[1]])))
u2[,1]=matrix(t(datos),D,1)

datos=mvrnorm(d,ceros,((phi2[2])*(a[[1]][[2]])))
u2[,2]=matrix(t(datos),D,1)
u11=matrix(0,D,k-1)
jj=1
for (i in 1:d){
	for(j in 1:t){
			u11[jj,]=u1[i,]
			jj=jj+1}}

x1=matrix(t(x1),d*t,1)
x2=matrix(t(x2),d*t,1,byrow=TRUE)
ind=matrix(rep(1.3,D),D,1)
ind2=matrix(rep(-1.6,D),D,1)
beta=c(-1,1)
pr=matrix(0,D,k-1)
theta=matrix(0,D,k-1)
for (j in 1:(k-1)){
	if (j==1) {theta[,j]=ind+x1*beta[j]+u11[,j]+u2[,j]}
	if (j==2) {theta[,j]=ind2+x2*beta[j]+u11[,j]+u2[,j]}
}

suma=rowSums(exp(theta))
a=1/(1+suma) 
for (i in 1:(k-1)){
	pr[,i]=a*exp(theta[,i])}
aa=list()
j=5
for ( i in 1:d){
	aa[[i]]=matrix(rep(j,t),t,1)
	j=j+5}
nu=do.call(rbind,aa)
aa=list()
j=200
for ( i in 1:d){
	aa[[i]]=matrix(rep(j,t),t,1)
	j=j+100}
nuu=do.call(rbind,aa)
y=matrix(0,D,(k))
yr=matrix(0,D,(k))
for (i in 1:D){
	y[i,]=t(rmultinom(1,nu[i],c(pr[i,1],pr[i,2],a[i])))
	yr[i,]=t(rmultinom(1,nuu[i]-nu[i],c(pr[i,1],pr[i,2],a[i])))}
a=list()
for ( i in 1:d){
	a[[i]]=matrix(rep(i,t),t,1)}
area=do.call(rbind,a)
time=rep(seq(1:t),d)
output=cbind(area,time,nu,nuu,y,cbind(x1,x2),yr)
return(output)}

#DATA
data=simulation(50,10,3)
colnames(data)=c("area","time","sample","population","y1","y2","y3","x1","x2","y11","y22","y33")
data=as.data.frame(data)
names(data)
data=subset(data,data$time>4)


library(mme)

k=3 #number of categories of the response variable
pp=c(1,1) #vector with the number of auxiliary variables in each category #data
mod=2 #Model 2
#Needed matrix and initial values
datar=data.mme(data[,1:9],k,pp,mod)


#Model fit
result=model(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],datar$n,datar$N, mod)
result

#Fixed effects
result$beta.Stddev.p.value

#Random effects
result$phi.Stddev.p.value


#Direct estimators
dir1=data$y11
dir2=data$y22


#Plot direct estimator in front of model estimator
dos.ver<-matrix(1:2,1,2)
layout(dos.ver)
plot(dir1,result$mean[,1],main="Small area estimator Y1",xlab="Direct estimate", ylab="model estimate",font.main=2,cex.main=1.5,cex.lab=1.3)
abline(a=0,b=1)
plot(dir2,result$mean[,2],main="Small area estimator Y2",xlab="Direct estimate", ylab="model estimate",font.main=2,cex.main=1.5,cex.lab=1.3)
abline(a=0,b=1)


#Model estimator
data$yest1=result$mean[,1]
data$yest2=result$mean[,2]


#Plot direct estimator and model estimator ordered by sample size for time=10
dos.ver<-matrix(1:2,1,2)
layout(dos.ver)

a=subset(data,data[,2]==10)
a=a[order(a[,3]),]
g_range <- range(0,45)
plot(a$y11/1000,type="b", col="blue",axes=FALSE, ann=FALSE)
lines(a$yest1/1000,type="b",pch=4, lty=2, col="red")
title(xlab="Sample size")
axis(1,at=c(1,10,20,30,40,50),lab=c(a$sample[1],a$sample[10],a$sample[20],a$sample[30],a$sample[40],a$sample[50]))
axis(2, las=1, at=1*0:g_range[2])
legend("topleft", c("Direct","Model"), cex=1, col=c("blue","red"), 
   lty=1:2,pch=c(1,4), bty="n")
title(main="Small area estimator Y1", font.main=1.2,cex.main=1)

plot(a$y22/1000,type="b",col="blue",axes=FALSE, ann=FALSE)
lines(a$yest2/1000,type="b",pch=4, lty=2, col="red")
title(xlab="Sample size")
axis(1,at=c(1,10,20,30,40,50),lab=c(a$sample[1],a$sample[10],a$sample[20],a$sample[30],a$sample[40],a$sample[50]))
axis(2, las=1, at=1*0:g_range[2])
legend("topleft", c("Direct","Model"), cex=1, col=c("blue","red"), 
   lty=1:2,pch=c(1,4), bty="n")
title(main="Small area estimator Y2", font.main=1.2,cex.main=1)

##Bootstrap parametric BIAS and MSE

B=100   #Bootstrap iterations
ss=12345 #SEED
set.seed(ss)

mse.pboot=mseb(pp,datar$Xk,datar$X,datar$Z,datar$n,datar$N,result,B,mod)

#RMSE
data$rmse1=mse.pboot[[3]][,1]
data$rmse2=mse.pboot[[3]][,2]


#PLOT THE RMSE ORDERED BY SAMPLE SIZE FOR TIME=10

a=subset(data,data[,2]==10)
a=a[order(a[,3]),]
dos.ver<-matrix(1:2,1,2)
layout(dos.ver)
g_range <- range(0,45)
plot(a$rmse1,type="b", col="blue",axes=FALSE, ann=FALSE)
title(xlab="Sample size")
axis(1,at=c(1,10,20,30,40,50),lab=c(a$sample[1],a$sample[10],a$sample[20],a$sample[30],a$sample[40],a$sample[50]))
axis(2, las=1, at=10*0:g_range[2])
title(main="RMSE for the estimator of Y1", font.main=1.2,cex.main=1)

g_range <- range(0,45)
plot(a$rmse2,type="b",col="blue",axes=FALSE, ann=FALSE)
title(xlab="Sample size")
axis(1,at=c(1,10,20,30,40,50),lab=c(a$sample[1],a$sample[10],a$sample[20],a$sample[30],a$sample[40],a$sample[50]))
axis(2, las=1, at=10*0:g_range[2])
title(main="RMSE for the estimator of Y2", font.main=1.2,cex.main=1)

