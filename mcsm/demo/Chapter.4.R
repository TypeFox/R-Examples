# Section 4.2, Monitoring variation

h=function(x){(cos(50*x)+sin(20*x))^2}
x=matrix(h(runif(200*10^4)),ncol=200)
estint=apply(x,2,cumsum)/(1:10^4)
plot(estint[,1],ty="l",col=0,ylim=c(.8,1.2),xlab="Iterations",ylab="")
y=apply(estint,1,quantile,c(.025,.975))
polygon(c(1:10^4,10^4:1),c(y[1,],rev(y[2,])),col="wheat",border=FALSE)

boot=matrix(sample(x[,1],200*10^4,rep=T),nrow=10^4,ncol=200)
bootit=apply(boot,2,cumsum)/(1:10^4)
bootup=apply(bootit,1,quantile,.975)
bootdo=apply(bootit,1,quantile,.025)
polygon(c(1:10^4,10^4:1),c(bootup,rev(bootdo)),col="gold3",border=FALSE)

lines(y[1,],lty=2,lwd=2)  	#bottom boundary

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Normal version

norma=matrix(rnorm(500*10^4),ncol=500)+2.5
weit=1/(1+norma^2)
esti=apply(norma*weit,2,cumsum)/apply(weit,2,cumsum)
plot(esti[,1],type="l",col="white",ylim=c(1.7,1.9),xlab="Iterations",ylab="")
band=apply(esti,1,quantile,c(.025,.975))
polygon(c(1:10^4,10^4:1),c(band[1,],rev(band[2,])),col="grey",border=FALSE)

vare=cumsum(weit[,1]*norma[,1]^2)/cumsum(weit[,1])-esti[,1]^2
lines(esti[,1]+2*sqrt(vare/(1:10^4)))
lines(esti[,1]-2*sqrt(vare/(1:10^4)))

varw=cumsum(weit[,1]^2)*(1:10^4)/cumsum(weit[,1])^2
lines(esti[,1]+2*sqrt(varw*vare/(1:10^4)),col="sienna")
lines(esti[,1]-2*sqrt(varw*vare/(1:10^4)),col="sienna")

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

#Cauchy version

cocha=matrix(rcauchy(500*10^4),ncol=500)
wach=dnorm(cocha,mean=2.5)
wachd=wach #Pb with zero weights
print(range(cocha))
print(range(wach))

wachd[wachd<10^(-6)]=10^(-6)
echi=apply(cocha*wach,2,cumsum)/apply(wachd,2,cumsum)
vare=cumsum(wach[,1]*cocha[,1]^2)/cumsum(wachd[,1])-echi[,1]^2
varw=cumsum(wachd[,1]^2)*(1:10^4)/cumsum(wachd[,1])^2
sek=seq(1,10^4,le=10^3)

par(mar=c(4,2,2,2))
plot(echi[,1],type="l",col="white",ylim=c(1.65,1.85),xlab="Iterations",ylab="")
polygon(c(1:10^4,10^4:1),c(band[1,],rev(band[2,])),bor=NA,col="grey")
lines(sek,(echi[,1]+2*sqrt(vare/(1:10^4)))[sek],lwd=2)
lines(sek,(echi[,1]-2*sqrt(vare/(1:10^4)))[sek],lwd=2)
lines(sek,(echi[,1]+2*sqrt(varw*vare/(1:10^4)))[sek],col="sienna",lwd=2)
lines(sek,(echi[,1]-2*sqrt(varw*vare/(1:10^4)))[sek],col="sienna",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 4.3, Effective sample size

ess=apply(weit,2,cumsum)^2/apply(weit^2,2,cumsum)
essbo=apply(ess,1,quantile,c(.025,.975))
ech=apply(wachd,2,cumsum)^2/apply(wachd^2,2,cumsum)
echbo=apply(ech,1,quantile,c(.025,.975))

par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(ess[,1],type="l",col=0,xlab="Iterations",ylab="ESS")
polygon(c(1:10^4,10^4:1),c(essbo[1,],rev(essbo[2,])),col="grey",bo=NA)
polygon(c(1:10^4,10^4:1),c(echbo[1,],rev(echbo[2,])),col="wheat",bo=NA)

sumweit=apply(weit,2,cumsum)
plex=(apply(weit*log(weit),2,cumsum)/sumweit)-log(sumweit)
chumweit=apply(wachd,2,cumsum)
plech=(apply(wachd*log(wachd),2,cumsum)/chumweit)-log(chumweit)
plob=apply(exp(plex),1,quantile,c(.025,.975))
ploch=apply(exp(plech),1,quantile,c(.025,.975))

plot(plex[,1]/(1:10^4),type="l",col=0,xlab="Iterations",ylab="Perplexity",ylim=c(.3,1))
polygon(c(1:10^4,10^4:1),c(plob[1,]/(1:10^4),rev(plob[2,]/(1:10^4))),col="grey",bo=NA)
polygon(c(1:10^4,10^4:1),c(ploch[1,]/(1:10^4),rev(ploch[2,]/(1:10^4))),col="wheat",bo=NA)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 4.4, Brownian bands

dev.off()
ustar=function(t){ 0.3+2.35*sqrt(t) }
curve(ustar,from=0,to=1,ylim=c(-3,3),ylab="",lwd=2)
curve(-ustar(x),from=0,to=1,add=T,lwd=2)
lines(c(0,0),c(-.3,.3),lwd=2)
lines(seq(0,1,le=10^4),cumsum(rnorm(10^4,sd=10^(-2))),col="steelblue",lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 4.4, Example 4.5

N=10^4;norma=rnorm(N)+2.5
hnorm=norma*dcauchy(norma)
munorm=mean(hnorm);sdnorm=sd(hnorm)

X11(h=3.5);par(mfrow=c(1,2),mar=c(4,2,2,1))
curve(munorm+(.1+3.15*sqrt(x))*sdnorm*10^2/round(x*N),
lwd=2,from=0,to=1,ylim=c(.06,.13),xlab="Iterations",ylab="")
curve(munorm-((.1+3.15*sqrt(x))*sdnorm*10^2/round(x*N)),
lwd=2,from=0,to=1,add=T)
f=function(x){(cumsum(hnorm))[round(N*x)]/round(x*N)}
curve(f,lwd=2,from=0.001,to=1,col="steelblue",add=T)

cocha=rcauchy(N)
hcoch=cocha*dnorm(cocha-2.5)
mucoch=mean(hcoch);sdcoch=sd(hcoch)

curve(mucoch+(.1+3.15*sqrt(x))*sdcoch*10^2/round(x*N),
lwd=2,from=0,to=1,add=T)
curve(mucoch-((.1+3.15*sqrt(x))*sdcoch*10^2/round(x*N)),
lwd=2,from=0,to=1,add=T)
f=function(x){(cumsum(hcoch))[round(N*x)]/round(x*N)}
curve(f,lwd=2,from=0.001,to=1,col="steelblue",add=T,lty=2)

hnorm=dcauchy(norma)
munorm=mean(hnorm);sdnorm=sd(hnorm)

curve(munorm+(.1+3.15*sqrt(x))*sdnorm*10^2/round(x*N),
lwd=2,from=0,to=1,ylim=c(.04,.1),xlab="Iterations",ylab="")
curve(munorm-((.1+3.15*sqrt(x))*sdnorm*10^2/round(x*N)),
lwd=2,from=0,to=1,add=T)
f=function(x){(cumsum(hnorm))[round(N*x)]/round(x*N)}
curve(f,lwd=2,from=0.001,to=1,col="steelblue",add=T)

hcoch=dnorm(cocha-2.5)
mucoch=mean(hcoch);sdcoch=sd(hcoch)

curve(mucoch+(.1+3.15*sqrt(x))*sdcoch*10^2/round(x*N),
lwd=2,from=0,to=1,add=T)
curve(mucoch-((.1+3.15*sqrt(x))*sdcoch*10^2/round(x*N)),
lwd=2,from=0,to=1,add=T)
f=function(x){(cumsum(hcoch))[round(N*x)]/round(x*N)}
curve(f,lwd=2,from=0.001,to=1,col="steelblue",add=T,lty=2)

S=readline(prompt="Type  <Return>   to continue : ")
dev.off();dev.off()

# Section 4.5, Rao-Blackwellization

nu=5;mu=3;sigma=0.5
y=sqrt(rchisq(N,df=nu)/nu)
x=rnorm(N,mu,sigma/y)
d1=cumsum(exp(-x^2))/(1:N)
d2=cumsum(exp(-mu^2/(1+2*(sigma/y)^2))/sqrt(1+2*(sigma/y)^2))/(1:N)
par(mar=c(4,2,2,1))
plot(d1,ty="l",lwd=2,xlab="Iterations",ylab="",ylim=c(.006,.009))
lines(d2,col="gold3",lwd=2)

S=readline(prompt="Type  <Return>   to continue (warning, lengthy step!): ")
dev.off()

# Section 4.6, Correlated simulations
jamestein()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 4.6, antithetic variables
dyadic()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 4.6, control variables

x=2.5
thet=rnorm(10^3,mean=x)
delt=thet/(1+thet^2)
moms=delta=c()
for (i in 1:5){
   moms=rbind(moms,(thet-x)^(2*i-1))
   reg=lm(delt~t(moms)-1)$coef
   delta=rbind(delta,as.vector(delt-reg%*%moms))
   }
plot(cumsum(delt)/(1:10^3),ty="l",xlab="Iterations",ylab="")
for (i in 1:5) lines(cumsum(delta[i,])/(1:10^3))

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()

# Section 4.6, Pima Indians logistic regression
logima()

S=readline(prompt="Type  <Return>   to continue : ")
dev.off()
