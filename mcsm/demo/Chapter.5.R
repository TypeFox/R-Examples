# Section 5.1, Numerical approximation

ref=rcauchy(5001)

f=function(y){-sum(log(1+(x-y)^2))}

mi=NULL
for (i in 1:400){
  x=ref[1:i]
  aut=optimise(f,interval=c(-10,10),max=T)
  mi=c(mi,aut$max)
  }

f=function(y){1/prod((1+(x-y)^2))}

mip=NULL
for (i in 1:400){
  x=ref[1:i]
  aut=optimise(f,interval=c(-10,10),max=T)
  mip=c(mip,aut$max)
  }

par(mfrow=c(1,2),mar=c(2,4,2,1))
plot(mi,ty="l",lwd=2,xlab="",ylab="arg")
lines(mip,col="sienna",lwd=2)

f=function(y){-sin(y*100)^2-sum(log(1+(x-y)^2))}

mig=NULL
for (i in 1:400){
  x=ref[1:i]
  aut=optimise(f,interval=c(-10,10),max=T)
  mig=c(mig,aut$max)
  }

f=function(y){exp(-sin(y*100)^2)/prod((1+(x-y)^2))}

mis=NULL

for (i in 1:401){
  x=ref[1:i]
  aut=optimise(f,interval=c(-10,10),max=T)
  mis=c(mis,aut$max)
  }
plot(mig,ty="l",lwd=2,xlab="",ylab="arg")
lines(mis,col="sienna",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# normal mixture in Example 5.2

#data
da=sample(rbind(rnorm(10^2),2.5+rnorm(3*10^2)))

#minus the log-likelihood function
like=function(mu){
  -sum(log((.25*dnorm(da-mu[1])+.75*dnorm(da-mu[2]))))
  }

#log-likelihood surface
mu1=mu2=seq(-2,5,le=250)
lli=matrix(0,nco=250,nro=250)

for (i in 1:250)
for (j in 1:250)
  lli[i,j]=like(c(mu1[i],mu2[j]))

par(mar=c(4,4,1,1))
image(mu1,mu2,-lli,xlab=expression(mu[1]),ylab=expression(mu[2]))
contour(mu1,mu2,-lli,nle=100,add=T)

#sequence of NR maxima
starts=matrix(c(1,1,-1,-1,4,-1,4,2.5,1,4,-1,4.5),nrow=2)

for (j in 1:dim(starts)[2]){
  sta=starts[,j]
  mmu=sta
  for (i in 1:(nlm(like,sta)$it))
    mmu=rbind(mmu,nlm(like,sta,iter=i)$est)
  lines(mmu[,1],mmu[,2],lwd=2)
  }

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.2.1, Stochastic search

h=function(x){(cos(50*x)+sin(20*x))^2}
rangom=h(matrix(runif(10^6),ncol=10^3))
monitor=t(apply(rangom,1,cummax))
plot(monitor[1,],type="l",col="white",xlab="Iterations",ylab=expression(h(theta)))
polygon(c(1:10^3,10^3:1),c(apply(monitor,2,max),rev(apply(monitor,2,min))),col="grey")
abline(h=optimise(h,int=c(0,1),max=T)$ob)

S=readline(prompt="Type  <Return>   to continue (warning, lengthy step!): ")

dev.off()
# Section 5.2.1, Cauchy example 5.4
ref=.1*rnorm(5001)
f=function(y){-y^2*.1-sum(log(1+(x-y)^2))}

maxv=loc=truv=trul=NULL

for (i in seq(10,length(ref),le=50)){

  prop=runif(10^3,-5,5)
  x=ref[1:i];
  tru=optimise(f,c(-5,5),max=T)
  trul=c(trul,tru$max)
  truv=c(truv,tru$ob)
  vale=apply(as.matrix(prop),1,f)
  loc=c(loc,prop[order(-vale)[1]])
  maxv=c(maxv,max(vale))
  }

par(mar=c(4,4,1,1),mfrow=c(2,1))
plot(trul,loc,cex=.5,pch=19,xlab=expression(theta^0),ylab=expression(hat(theta)))
abline(a=0,b=1,col="grey")
plot(seq(10,length(ref),le=50),(truv-maxv)/abs(truv),type="l",lwd=2,xlab="Sample size",ylab="Relative error")

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.2.1, Cauchy example 5.5
N=1+10^2
dft=(N-1)/2
cau=rcauchy(N)
mcau=median(cau)
rcau=diff(quantile(cau,c(.25,.75)))
f=function(x){
  z=dcauchy(outer(x,cau,FUN="-"));apply(z,1,mean)}
fcst=integrate(f,low=-20,up=20)$va
g=function(x){dt((x-mcau)/rcau,df=dft)/rcau}
ft=function(x){f(x)/fcst}

par(mar=c(2,2,1,1))
curve(ft,from=-10,to=10,xlab="",ylab="",lwd=2,n=10^3);curve(g,add=T,col="steelblue",lwd=2,lty=2,n=10^3)

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.2.1, Cauchy example 5.5 cont.
unisan=matrix(f(runif(500*10^2,-5,5)),ncol=500)
causan=matrix(f(rt(500*10^2,df=dft)*rcau+mcau),ncol=500)
unimax=apply(unisan,2,cummax)[10:10^2,]
caumax=apply(causan,2,cummax)[10:10^2,]
plot(caumax[,1],col="white",ylim=c(.8*max(causan),max(causan)))
polygon(c(10:10^2,10^2:10),c(apply(unimax,1,max),rev(apply(unimax,1,min))),col="grey")
polygon(c(10:10^2,10^2:10),c(apply(caumax,1,max),rev(apply(caumax,1,min))),col="wheat")


S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.2.1, artificial example 5.6:
#highly multimodal function to maximise
h=function(x){-(x[1]*sin(20*x[2])+x[2]*sin(20*x[1]))^2*cosh(sin(10*x[1])*x[1])-
 (x[1]*cos(10*x[2])-x[2]*sin(10*x[1]))^2*cosh(cos(20*x[2])*x[2])}

ha=function(x,y){-(x*sin(20*y)+y*sin(20*x))^2*cosh(sin(10*x)*x)-
(x*cos(10*y)-y*sin(10*x))^2*cosh(cos(20*y)*y)}

x=y=seq(-3,3,le=435)
z=-outer(x,y,ha)
par(bg = "wheat",mar=c(1,1,1,1))
persp(x, y, z, theta=155, phi=30, col="green4",
  ltheta=-120, shade=0.75, border =NA, box=FALSE)

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.2.2, artificial example 5.7
#Stochatic gradient version
start=c(.65,.8)
theta=matrix(start,ncol=2)
dif=iter=alpha=beta=1
hcur=hval=h(start)

while (dif>10^-5){

  zeta=rnorm(2);zeta=zeta/sqrt(t(zeta)%*%zeta)
  grad=alpha*zeta*(h(theta[iter,]+beta*zeta)-
        h(theta[iter,]-beta*zeta))/beta
  #safety condition  
  scale=sqrt(t(grad)%*%grad)
  while (scale>1){
    zeta=rnorm(2);zeta=zeta/sqrt(t(zeta)%*%zeta)
    grad=alpha*zeta*(h(theta[iter,]+beta*zeta)-
          h(theta[iter,]-beta*zeta))/beta
    scale=sqrt(t(grad)%*%grad)}
  theta=rbind(theta,theta[iter,]+grad)
  dif=scale;iter=iter+1
  alpha=1/(iter+1);beta=1/sqrt(iter+1)
  }

x=y=seq(-1,1,le=435)
z=outer(x,y,ha)
image(x,y,z,col=terrain.colors(150))
lines(theta,lwd=2)
points(theta[1,1],theta[1,2],col="gold",pch=19)
title(main=paste("min",format(-max(hval),dig=3),sep=" "))

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.2.2, artificial function example 5.8

h=function(x){(cos(50*x)+sin(20*x))^2}
par(mar=c(4,4,1,1),mfrow=c(2,2))
for (tt in 1:4){

  curve(h,from=0,to=1,n=10001,col="grey",lwd=2)
  sam=maximple()
  xgo=sam$x
  hgo=sam$y
  lines(xgo,hgo,col="steelblue4",lwd=2)
  points(xgo,hgo,col="steelblue4",cex=.5,pch=19)
  }

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.2.2, normal mixture example 5.9

da=sample(rbind(rnorm(10^2),2.5+rnorm(3*10^2)))

#minus the log-likelihood function
like=function(mu){
  -sum(log((.25*dnorm(da-mu[1])+.75*dnorm(da-mu[2]))))
  }

#log-likelihood surface
mu1=mu2=seq(-2,5,le=250)
lli=matrix(0,nco=250,nro=250)

for (i in 1:250)
for (j in 1:250)
  lli[i,j]=like(c(mu1[i],mu2[j]))


par(mar=c(4,4,1,1))
image(mu1,mu2,-lli,xlab=expression(mu[1]),ylab=expression(mu[2]))
contour(mu1,mu2,-lli,nle=100,add=T)

starts=matrix(c(1,1,-1,-1,4,-1,4,2.5,1,4,-1,4.5),nrow=2)

for (j in 1:dim(starts)[2]){
    sar=SAmix(starts[,j])
    lines(sar$the[,1],sar$the[,2],lwd=2,col="tomato");points(sar$the[sar$it,1],sar$the[sar$it,2],cex=.5,pch=19)
    }

# Section 5.2.3, artificial example 5.10
#SA version
start=c(.65,.8)
the=matrix(start,ncol=2)
temp=scale=dif=iter=factor=1
hcur=hval=h(start)

while (dif>10^-5){

  prop=the[iter,]+scale[iter]*rnorm(2)
  while (max(abs(prop))>1) prop=the[iter,]+scale[iter]*rnorm(2)

  if (temp[iter]*log(runif(1))>h(prop)-hcur)
    prop=the[iter,]

  hcur=h(prop);the=rbind(the,prop);hval=c(hval,hcur)
  iter=iter+1;temp=c(temp,1/(1*log(1+iter)));
  ace=length(unique(the[(iter/2):iter,1]))
  if (iter>100){
    if (ace==1) factor=factor/10
    if (2.5*ace>iter) factor=min(.1,factor*10)
    }
  scale=c(scale,min(.1,5*factor*sqrt(temp[iter])))
  dif=(iter<10^3)+(ace<2)+(max(hval)-max(hval[1:(iter/2)]))
  }

x=y=seq(-1,1,le=435)
z=outer(x,y,ha)
image(x,y,z,col=terrain.colors(150))
lines(the,lwd=2)
points(prop[1],prop[2],col="gold",pch=19)
title(main=paste("min",format(-max(hval),dig=3),sep=" "))

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.3.1, Pima indian example 5.11
pimax(10^2)

S=readline(prompt="Type  <Return>   to continue (Warning, lengthy run!): ")

dev.off()
# Section 5.3.3, EM algorithm
EMcenso()

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
# Section 5.3.4, MCEM algorithm
randogit(10^3,10^2)

S=readline(prompt="Type  <Return>   to continue (Warning, lengthy run!): ")
dev.off()
