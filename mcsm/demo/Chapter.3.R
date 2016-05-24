# Section 3.1, Introduction

ch=function(la){ integrate(function(x){x^(la-1)*exp(-x)},0,Inf)$val}
plot(lgamma(seq(.01,10,le=100)),log(apply(as.matrix(
 seq(.01,10,le=100)),1,ch)),xlab="log(integrate(f))",
 ylab=expression(log(Gamma(lambda))),pch=19,cex=.6)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

cac=rcauchy(10)+350
lik=function(the){u=dcauchy(cac[1]-the);for (i in 2:10)
 u=u*dcauchy(cac[i]-the);return(u)}

print(integrate(lik,-Inf,Inf))
print(integrate(lik,200,400))

library(MASS)
cac=rcauchy(10)
nin=function(a){integrate(lik,-a,a)$val}
nan=function(a){area(lik,-a,a)}
x=seq(1,10^3,le=10^4)
y=log(apply(as.matrix(x),1,nin))
z=log(apply(as.matrix(x),1,nan))
plot(x,y,type="l",ylim=range(cbind(y,z)),lwd=2)
lines(x,z,lty=2,col="sienna",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 3.2, Monte Carlo integration

h=function(x){(cos(50*x)+sin(20*x))^2}
par(mar=c(2,2,2,1),mfrow=c(2,1))
curve(h,xlab="Function",ylab="",lwd=2)
print(integrate(h,0,1))

x=h(runif(10^4))
estint=cumsum(x)/(1:10^4)
esterr=sqrt(cumsum((x-estint)^2))/(1:10^4)
plot(estint, xlab="Mean and error range",type="l",lwd=
 	2,ylim=mean(x)+20*c(-esterr[10^4],esterr[10^4]),ylab="")
lines(estint+2*esterr,col="gold",lwd=2)
lines(estint-2*esterr,col="gold",lwd=2)

S=readline(prompt="Type  <Return>   to continue (warning, lenghty step!) : ")

dev.off()
x=rnorm(10^7) #whole sample
bound=qnorm(c(.5,.75,.8,.9,.95,.99,.999,.9999))
res=matrix(0,ncol=8,nrow=8)
for (i in 1:8) #lengthy loop!!
for (j in 1:8)
   res[i,j]=mean(x[1:10^(i-1)]<bound[j])
print(matrix(as.numeric(format(res,digi=4)),ncol=8))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 3.3, Importance sampling

m=10^3
y=rexp(m)+4.5
weit=dnorm(y)/dexp(y-4.5)
plot(cumsum(weit)/1:m,type="l",xlab="Iterations",ylab="",lwd=2)
abline(a=pnorm(-4.5),b=0,col="gold3",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 3.3, Importance sampling, Example 3.6

f=function(a,b){
    exp(2*(lgamma(a+b)-lgamma(a)-lgamma(b))+a*log(.3)+b*log(.2))
    }
aa=1:150
bb=1:100
post=matrix(0,ncol=100,nrow=150)
for (i in 1:150)
  for (j in 1:100)
    post[i,j]=f(aa[i],bb[j])

X11(h=3.5);par(mfrow=c(1,2),mar=c(4,4,2,1))
image(aa,bb,post,xlab=expression(alpha),ylab=expression(beta))
contour(aa,bb,post,add=T)

x=matrix(rt(2*10^4,3),ncol=2)       #T sample
E=matrix(c(220,190,190,180),ncol=2) #Scale matrix
image(aa,bb,post,xlab=expression(alpha),ylab=expression(beta))
y=t(t(chol(E))%*%t(x)+c(50,45))
points(y,cex=.6,pch=19)

ine=apply(y,1,min);y=y[ine>0,];x=x[ine>0,]
normx=sqrt(x[,1]^2+x[,2]^2)
f=function(a) exp(2*(lgamma(sum(a))-sum(lgamma(a)))+ a[1]*log(.3)+a[2]*log(.2))
h=function(a) exp(1*(lgamma(sum(a))-sum(lgamma(a)))+ a[1]*log(.5)+a[2]*log(.5))
den=apply(as.matrix(normx),1,dt,3)
mean(apply(y,1,f)/den)/mean(apply(y,1,h)/den)
mean(y[,1]*apply(y,1,f)/den)/mean(apply(y,1,h)/den)
mean(y[,2]*apply(y,1,f)/den)/mean(apply(y,1,h)/den)

S=readline(prompt="Type  <Return>   to continue (warning, lengthy step!): ")
dev.off()

# Section 3.3.2, Sampling importance resampling (Example 3.6 continued)

par(mfrow=c(2,2),mar=c(4,4,2,1))
weit=(apply(y,1,f)/den)/mean(apply(y,1,h)/den)
image(aa,bb,post,xlab=expression(alpha),ylab=expression(beta))
points(y[sample(1:length(weit),10^3,rep=T,pro=weit),], cex=.6,pch=19)
boxplot(weit,ylab="importance weight")
plot(cumsum(weit)/(1:length(weit)),type="l",xlab="simulations", ylab="marginal likelihood")
boot=matrix(0,ncol=length(weit),nrow=100)
for (t in 1:100) boot[t,]=cumsum(sample(weit))/(1:length(weit))
uppa=apply(boot,2,quantile,.95);lowa=apply(boot,2,quantile,.05)
polygon(c(1:length(weit),length(weit):1),c(uppa,rev(lowa)), col="gold")
lines(cumsum(weit)/(1:length(weit)),lwd=2)
plot(cumsum(weit)^2/cumsum(weit^2),type="l", xlab="simulations", ylab="Effective sample size",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 3.3.3, Selection of the importance function

x=rnorm(10^6)
wein=dcauchy(x)/dnorm(x)
X11(h=3.5);par(mfrow=c(1,2),mar=c(4,4,2,1))
boxplot(wein/sum(wein),ylab="importance weight")
plot(cumsum(wein*(x>2)*(x<6))/cumsum(wein),type="l",xlab="Iterations",ylab="P(2<X<6)",lwd=2)
abline(a=pcauchy(6)-pcauchy(2),b=0,col="sienna",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 3.3.3, Selection of the importance function

h=function(x){(x>1)/sqrt(abs(x-1))}
sam1=rt(.95*10^4,df=2)
sam2=1+.5*rt(.05*10^4,df=2)^2
sam=sample(c(sam1,sam2),.95*10^4)
weit=dt(sam,df=2)/(0.95*dt(sam,df=2)+.05*(sam>0)*dt(sqrt(2*abs(sam-1)),df=2)*sqrt(2)/sqrt(abs(sam-1)))
plot(cumsum(h(sam1))/(1:length(sam1)),ty="l",xlab="Iterations",ylab="",lwd=2)
lines(cumsum(weit*h(sam))/1:length(sam1),col="blue",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 3.3.3, Probit modelling of Pima Indians

glm(type~bmi,data=Pima.tr,family=binomial(link="probit"))

post=function(be1,be2){
 mia=mean(Pima.tr$bmi)
 prod(pnorm(be1+(Pima.tr$bm[Pima.tr$t=="Yes"]-mia)*be2))*
 prod(pnorm(-be1-(Pima.tr$bm[Pima.tr$t=="No"]-mia)*be2))/exp((be1^2+be2^2)/200)
 }

ppost=function(be) post(be[1],be[2])

be1=seq(-.6,-.3,le=120);be2=seq(.04,.09,le=135)
lval=matrix(0,nrow=length(be1),ncol=length(be2))
for (i in 1:length(be1)) for (j in 1:length(be2)) lval[i,j]=post(be1[i],be2[j])
image(be1,be2,lval,xlab=expression(beta[1]),ylab=expression(beta[2]))

like=function(beda){
  mia=mean(Pima.tr$bmi)
  prod(pnorm(beda[1]+(Pima.tr$bm[Pima.tr$t=="Yes"]-mia)*beda[2]))*
  prod(pnorm(-beda[1]-(Pima.tr$bm[Pima.tr$t=="No"]-mia)*beda[2]))/
  exp(sum(beda^2)/200)}

sim=cbind(rnorm(10^3,m=-.4,sd=.04),rnorm(10^3,m=0.065,sd=.005))
weit=apply(sim,1,like)/(dnorm(sim[,1],m=-.4,sd=.04)*dnorm(sim[,2],m=0.065,sd=.005))
points(sim[sample(1:10^3,10^3,weit,rep=T),],pch=19,cex=.4)

sim=rbind(sim[1:(.95*10^3),],cbind(rnorm(.05*10^3,sd=10),rnorm(.05*10^3,sd=10)))
weit=apply(sim,1,ppost)/(.95*dnorm(sim[,1],m=-.4,sd=.081)*
 dnorm(sim[,2],m=0.065,sd=.01)+.05*dnorm(sim[,1],sd=10)*dnorm(sim[,2],sd=10))
X11();image(be1,be2,lval,xlab=expression(beta[1]),ylab=expression(beta[2]))
points(sim[sample(1:10^3,10^3,weit,rep=T),],pch=19,cex=.4)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off();dev.off()
