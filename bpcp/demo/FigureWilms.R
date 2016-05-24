
CItypes<-(c("Greenwood","Modified Lower","Beta Product","Borkowf","Strawderman-Wells","Thomas-Grunkemeier"))
CILtypes<-(c(1,1,1,1,4,3))
CILWD<-(c(3,1,1,1,1,1.5))
CIgray<-(c(.7,.7,0,.7,.7,.7))

library(bpcp)

## default oma=c(0,0,0,0)
par(oma=c(1,0,0,0))
library(survival)
## data from survival package
time<-nwtco$edrel/365.25
status<-nwtco$rel
d<-data.frame(time=time,status=status)
gw<-survfit(Surv(time,status)~1,data=d)
bfit<-bpcp(time,status)
plot(bfit,mark.time=FALSE,ylim=c(.7,1),lwd=1,linetype="surv",xlab="",axes=FALSE,ylab="Proportion Relapse-Free")
axis(1,at=2*0:8,label=2*0:8)
axis(2)
box(bty="l")
natrisk<-function(time,status,at=2*0:8){
    n<-rep(NA,length(at))
    for (i in 1:length(at)){
        n[i]<- length(time[time>at[i]])
    }
    list(at=at,n=n)
}
r<-natrisk(time,status)
mtext("Years",side=1,line=2,at=8)
mtext("Number at Risk",side=1,line=4.5,at=8,cex=.8)
mtext(r$n,side=1,at=r$at,line=3.5,cex=.8)
class(gw)<-"kmci"
lines(gw,col=gray(CIgray[1]),conf.int="only",lwd=CILWD[1],lty=CILtypes[1])
#lines(mod,col=gray(CIgray[2]),conf.int="only",lwd=CILWD[2],lty=CILtypes[2])
lines(bfit,col=gray(CIgray[3]),lwd=CILWD[3],lty=CILtypes[3])
#lines(borkowf,col=gray(CIgray[4]),lwd=CILWD[4],lty=CILtypes[4])
#lines(SW,col=gray(CIgray[5]),lwd=CILWD[5],lty=CILtypes[5])
#lines(TG,col=gray(CIgray[6]),lwd=CILWD[6],lty=CILtypes[6])

o<-c(3,1,2,4,5,6)
o<-c(3,1)
legend(1,.8,
   lty=CILtypes[o],
   lwd=CILWD[o],
   col=gray(CIgray[o]),
   CItypes[o],
   bty="n" )


