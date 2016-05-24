## sclerosis data set is in bpcp R package, 
## Wilms Tumor data set is in survival R package
library(bpcp)
library(survival)

##############################
#
#   Sclerosis Data 
###############################
par(mfrow=c(1,2),mar=c(5,4,4,2)+.1,oma=c(0,0,0,0))
## Sclerosis Data
## data from bpcp package
library(bpcp)
data(sclerosis)
time<-sclerosis$year
status<-sclerosis$status
bout<-bpcp(time,status,midp=FALSE)
plot(bout,ciLTY=1,ciCOL="gray",lwd=6,xlab="Year",ylab="Survival",ylim=c(0,1),mark.time=FALSE)
bout2<-bpcp(time,status,midp=FALSE,monotonic=FALSE)
lines(bout2,col="red",lwd=3,lty=1)
bout3<-bpcp(time,status,midp=FALSE,monotonic=FALSE,nmc=10^4)
lines(bout3,col="blue",lwd=1,lty=1)

legend("bottomleft",legend=c("MM (mono)","MM","MC (m=10,000)"),
    lwd=c(6,3,1),lty=c(1,1,1),col=c("gray","red","blue"))

title("Standard BPCP")

bmid<-bpcp(time,status,midp=TRUE)
plot(bmid,ciLTY=1,ciCOL="gray",lwd=6,xlab="Year",ylab="Survival",ylim=c(0,1),mark.time=FALSE)
#bmid<-bpcp(time,status,midp=TRUE)
#lines(bmid,col="black",lwd=1,lty=1)
bmid2<-bpcp(time,status,midp=TRUE,monotonic=FALSE)
lines(bmid2,col="red",lwd=3,lty=1)
bmid3<-bpcp(time,status,midp=TRUE,monotonic=FALSE,nmc=10^4)
lines(bmid3,col="blue",lwd=1,lty=1)

title("mid-p BPCP")

legend("bottomleft",legend=c("MM (mono)","MM","MC (m=10,000)"),
    lwd=c(6,3,1),lty=c(1,1,1),col=c("gray","red","blue"))


#dev.print(postscript,file="../tex/SSclerosis.eps",horizontal=FALSE)

##########################################
#
# MM vs MC plots
#   Sclerosis data 
##########################################



data(sclerosis)
sort(sclerosis$year[sclerosis$status==1])
time<-sclerosis$day
year<-sclerosis$year
Tyr<-sort(year[status==1])

status<-sclerosis$status
T<-sort(time[status==1])
Y<-function(t){
   # get number at risk just before t
   length(time[time>=t])
}

k<-length(T)
y<-rep(NA,k)
for (i in 1:k){
    y[i]<-Y(T[i])    
}
y

b<-bpcp(time,status,midp=FALSE,monotonic=FALSE)
b

aT<-bT<-rep(NA,k)
for (i in 1:k){
    aT[i]<-b$betaParms$alower[b$R==T[i] & b$L==T[i]]
    bT[i]<-b$betaParms$blower[b$R==T[i] & b$L==T[i]]
}

aT
bT

#abmm(aT[10],bT[10],y[11],1)



set.seed(1031)
simBB<-function(a1,b1,a2,b2,nsim=10^4,xtext=c(.2,.2),ytext=c(8,6),...){
    B1<-rbeta(nsim,a1,b1)
    B2<-rbeta(nsim,a2,b2)
    BB<-B1*B2
    d<-density(BB)
    mm<-abmm(a1,b1,a2,b2)
    p<-0:1000/1000
    dmm<-dbeta(p,mm$a,mm$b)
    par(mfrow=c(2,1))
    dB1<-dbeta(p,a1,b1)
    dB2<-dbeta(p,a2,b2)
    plot(c(p,p),c(dB1,dB2),type="n",xlab="",ylab="")
    lines(p,dB1,col="red",lwd=3)
    lines(p,dB2,col="blue",lwd=3)
    text(xtext[1],ytext[1],paste0("B(",round(a1,4),",",round(b1,4),")"),col="red")
    text(xtext[2],ytext[2],paste0("B(",round(a2,4),",",round(b2,4),")"),col="blue")


    plot(0,0,type="n",ylim=c(0,max(c(d$y,dmm))),xlab="",ylab="",...)
    lines(d$x,d$y,col=gray(.8),lwd=3)
    lines(p,dmm,lwd=3,lty=2)
    legend("topleft",legend=c("method of moments","Monte Carlo"),lwd=c(3,3),lty=c(2,1),col=gray(c(0,.8)))
    ## set back to default
    par(mfrow=c(1,1))

    return(d)
}
simBB(aT[10],bT[10],y[11],1,xlim=c(0,1))
Tyr[11]
T[11]
Tyr[12]
T[12]
#dev.print(postscript,file="../tex/MMvsMC1.eps",horizontal=FALSE)

simBB(aT[11],bT[11],y[12],1,xlim=c(0,1),ytext=c(3,2))

#dev.print(postscript,file="../tex/MMvsMC2.eps",horizontal=FALSE)

##############################
#
#   Wilms Tumor Data 
###############################
par(mfrow=c(1,2),mar=c(5,4,4,2)+.1,oma=c(0,0,0,0))

## NWTCO Data
## data from survival package
library(survival)
time<-nwtco$edrel/365.25
status<-nwtco$rel
bout<-bpcp(time,status,midp=FALSE)
plot(bout,ciLTY=1,ciCOL="gray",lwd=6,xlab="Year",
     ylab="Proportion Relapse-Free",ylim=c(0.7,1),mark.time=FALSE)
bout2<-bpcp(time,status,midp=FALSE,monotonic=FALSE)
lines(bout2,col="red",lwd=3,lty=1)
bout3<-bpcp(time,status,midp=FALSE,monotonic=FALSE,nmc=10^4)
lines(bout3,col="blue",lwd=1,lty=1)

legend("bottomleft",legend=c("MM (mono)","MM","MC (m=10,000)"),
    lwd=c(6,3,1),lty=c(1,1,1),col=c("gray","red","blue"))

title("Standard BPCP")

bmid<-bpcp(time,status,midp=TRUE)
plot(bmid,ciLTY=1,ciCOL="gray",lwd=6,xlab="Year",
    ylab="Proportion Relapse-Free",ylim=c(0.7,1),mark.time=FALSE)
#bmid<-bpcp(time,status,midp=TRUE)
#lines(bmid,col="black",lwd=1,lty=1)
bmid2<-bpcp(time,status,midp=TRUE,monotonic=FALSE)
lines(bmid2,col="red",lwd=3,lty=1)
bmid3<-bpcp(time,status,midp=TRUE,monotonic=FALSE,nmc=10^4)
lines(bmid3,col="blue",lwd=1,lty=1)

title("mid-p BPCP")

legend("bottomleft",legend=c("MM (mono)","MM","MC (m=10,000)"),
    lwd=c(6,3,1),lty=c(1,1,1),col=c("gray","red","blue"))


#dev.print(postscript,file="../tex/SWilms.eps",horizontal=FALSE)

##########################################
#
# MM vs MC plots
#   Wilms Tumor data 
##########################################
time<-nwtco$edrel
#/365.25
status<-nwtco$rel
T<-sort(unique(time[status==1]))
k<-length(T)
y<-aT<-bT<-rep(NA,k)
b<-bpcp(time,status,midp=FALSE,monotonic=FALSE)
for (i in 1:k){
    y[i]<-Y(T[i])   
    aT[i]<-b$betaParms$alower[b$R==T[i] & b$L==T[i]]
    bT[i]<-b$betaParms$blower[b$R==T[i] & b$L==T[i]] 
}

allt<-sort(unique(time))
length(allt)
bmidp<-bpcp(time,status,midp=TRUE,monotonic=FALSE)

picki<-function(i){
     I<-b$R==allt[i] & b$L==allt[i]
     list(Interval=b$Interval[I],
          t=b$L[I],
          alo=b$betaParms$alower[I],
          blo=b$betaParms$blower[I],
          ahi=b$betaParms$aupper[I],
          bhi=b$betaParms$bupper[I],
          umid=bmidp$upper[I],
          y=Y(b$L[I]))
}

z<-picki(2489)
z$t
T[k]
y[k]
Y(z$t)

midpBetaPlots<-function(alo,blo,ahi,bhi,umid,...){
     p<-0:1000/1000
     Blo<-dbeta(p,alo,blo)
     Bhi<-dbeta(p,ahi,bhi)
     plot(c(0,1),range(c(0,Blo,Bhi)),type="n",ylab="",xlab="",...)
     lines(p,Blo,col="gray",lwd=6)
     lines(p,Bhi,col="black",lwd=2,lty=1)
     ulo<-qbeta(.975,alo,blo)
     uhi<-qbeta(.975,ahi,bhi)
     lines(c(ulo,ulo),c(0,-5),col="gray",lwd=6)
     lines(c(uhi,uhi),c(0,-5),col="black",lwd=2,lty=1)
     lines(c(umid,umid),c(0,-5),col="red",lwd=2)
     return(list(ulo=ulo,uhi=uhi))
}
par(mfrow=c(1,1))
midpBetaPlots(z$alo,z$blo,z$ahi,z$bhi,z$umid,xlim=c(.81,.88))
legend("topleft",legend=c("lower BPCP Conf Distn (MM)","upper BPCP Conf Distn (MM)",
                      "upper 95% mid-p BPCP (MM)"),col=c("gray","black","red"),
                       lwd=c(6,2,2),lty=c(1,1,1))

#dev.print(postscript,file="../tex/midpPlot1.eps",horizontal=FALSE)

z<-picki(2730)
par(mfrow=c(1,1))
midpBetaPlots(z$alo,z$blo,z$ahi,z$bhi,z$umid,xlim=c(.7,.95))
legend("topleft",legend=c("lower BPCP Conf Distn (MM)","upper BPCP Conf Distn (MM)",
                      "upper 95% mid-p BPCP (MM)"),col=c("gray","black","red"),
                       lwd=c(6,2,2),lty=c(1,1,1))

#dev.print(postscript,file="../tex/midpPlot2.eps",horizontal=FALSE)


####################################
#
# Simulated Data
####################################
createDataBeta<-function(n,param.a,param.b){
  x<-rbeta(n,param.a,param.b)
  cens<-runif(n)
  time<-pmin(x,cens)
  status<-rep(0,n)
  status[time==x]<-1
  list(time=time,status=status,xi=x,ci=cens)
}
## checked seeds from 1:18, 18 is first to give non-monotonic upper limit
set.seed(18)
d<-createDataBeta(1000,.1,.1)
b<-bpcp(d$time,d$status,midp=FALSE,monotonic=FALSE)
b2<-bpcp(d$time,d$status,midp=FALSE,monotonic=TRUE)
par(mfrow=c(1,2))
plot(b2$R,b2$upper,type="l",lwd=6,col=gray(.8),
    xlab="time",ylab="Upper Limit from 95% BPCP",xlim=c(0,1),ylim=c(0,1),main="Simulated Data")
lines(b$R,b$upper,lwd=2,col="red")

legend("bottomleft",legend=c("MM,mono","MM"),
                col=c("gray","red"),
                       lwd=c(6,2),lty=c(1,1))


plot(b2$R,b2$upper,type="l",lwd=6,col=gray(.8),
     xlab="time",ylab="Upper Limit from 95% BPCP",xlim=c(.94,1),ylim=c(.4,.5),main="Simulated Data")
lines(b$R,b$upper,lwd=2,col="red")


legend("bottomleft",legend=c("MM,mono","MM"),
                col=c("gray","red"),
                       lwd=c(6,2),lty=c(1,1))


#dev.print(postscript,file="../tex/BPCPnonMono1.eps",horizontal=FALSE)

# get number at risk and failures
#  totalFailures
sum(d$status)
## 515th failure, just after tt1, where
tt1<-.9480465
b$R[b$R>tt1]
b$upper[b$R>tt1]
b$betaParms$aupper[b$R>tt1]
b$betaParms$bupper[b$R>tt1]


k<-kmgw.calc(d$time,d$status)
k$time[k$time>tt1]
k$ni[k$time>tt1]
k$di[k$time>tt1]
# last failure just after tt2
tt2<-0.9976
k$time[k$time>tt2]
k$ni[k$time>tt2]
k$di[k$time>tt2]

library(survival)
plot(survfit(Surv(d$time,d$status)~1))

a1<-b$betaParms$aupper[b$R>tt1][3]
b1<-b$betaParms$bupper[b$R>tt1][3]
a2<-b$betaParms$aupper[b$R>tt1][40]
b2<-b$betaParms$bupper[b$R>tt1][40]
# check that a2,b2 is the MM of B(a1,b1)B(2,1)
mm<-abmm(a1,b1,2,1)

p<-0:1000/1000
par(mfrow=c(1,1))
plot(p,dbeta(p,a1,b1),type="n",xlab="x",ylab="f(x)",main="Beta densities")
lines(p,dbeta(p,2,1),col="red",lwd=2)
lines(p,dbeta(p,a1,b1),col="gray",lwd=6)
lines(p,dbeta(p,a2,b2),col="blue",lwd=2)

set.seed(1)
simBBline<-function(a1,b1,a2,b2,nsim=10^5){
    B1<-rbeta(nsim,a1,b1)
    B2<-rbeta(nsim,a2,b2)
    BB<-B1*B2
    d<-density(BB)
    lines(d$x,d$y,col="green3",lwd=2)
    return(BB)
}

BB<-simBBline(a1,b1,2,1)



lines(rep(qbeta(.975,a1,b1),2),c(-1,-.1),col="gray",lwd=6)
lines(rep(qbeta(.975,a2,b2),2),c(-1,-.1),col="blue",lwd=2)
lines(rep(quantile(BB,probs=.975),2),c(-1,-.1),col="green3",lwd=2)

qbeta(.975,a1,b1)
qbeta(.975,a2,b2)
quantile(BB,probs=.975)
dev.print(postscript,file="../tex/BPCPnonMono2.eps",horizontal=FALSE)

