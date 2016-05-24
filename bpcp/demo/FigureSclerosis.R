

### see http://jfly.iam.u-toyo.ac.jp/color/#assign
colblind8<-rep(NA,8)
names(colblind8)<-c("black","orange","skyBlue","bluishGreen",
     "yellow","blue","vermillion","reddishPurple")
colblind8[1]<-rgb(0,0,0,maxColorValue=255)
## orange
colblind8[2]<-rgb(230,159,0,maxColorValue=255)
## sky blue
colblind8[3]<-rgb(86,180,233,maxColorValue=255)
### bluish green
colblind8[4]<-rgb(0,158,115,maxColorValue=255)
colblind8[5]<-rgb(240,228,66,maxColorValue=255)
colblind8[6]<-rgb(0,114,178,maxColorValue=255)
colblind8[7]<-rgb(213,94,0,maxColorValue=255)
colblind8[8]<-rgb(204,121,167,maxColorValue=255)


CItypes<-(c("Greenwood","Modified Lower","Beta Product","Borkowf","Strawderman-Wells","Thomas-Grunkemeier"))
CILtypes<-(c(2,2,1,1,4,3))
CILWD<-(c(2,1,2,1,1,1.5))
CIgray<-(c(.7,.7,0,.7,.7,.7))
CIcol<-colblind8[c(3,2,1,4,5,6)]
CIcol<-c(rep(gray(.7),2),gray(0),rep(gray(.7),3))

library(bpcp)
data(sclerosis)
library(survival)

mod<-survfit(Surv(year,status)~1,data=sclerosis,conf.lower='modified')
gw<-survfit(Surv(year,status)~1,data=sclerosis)
new<-bpcp(sclerosis$year,sclerosis$status)
borkowf<-kmciBorkowf(sclerosis$year,sclerosis$status)
SW<-kmciSW(sclerosis$year,sclerosis$status)
TG<-kmciTG(sclerosis$year,sclerosis$status)
plot(gw,lwd=2,xlab="years",ylab="survival",col=gray(1),las=1,xaxt="none",bty="l")
axis(1,at=0:8)
lines(gw,lwd=2)
lines(gw,col=CIcol[1],conf.int="only",lwd=CILWD[1],lty=CILtypes[1])
lines(mod,col=CIcol[2],conf.int="only",lwd=CILWD[2],lty=CILtypes[2])
lines(new,col=CIcol[3],lwd=CILWD[3],lty=CILtypes[3])
lines(borkowf,col=CIcol[4],lwd=CILWD[4],lty=CILtypes[4])
lines(SW,col=CIcol[5],lwd=CILWD[5],lty=CILtypes[5])
lines(TG,col=CIcol[6],lwd=CILWD[6],lty=CILtypes[6])

o<-c(3,1,2,4,5,6)
#o<-c(3,1)
legend(1,.4,
   lty=CILtypes[o],
   lwd=CILWD[o],
   col=CIcol[o],
   CItypes[o],
   bty="n" )



  