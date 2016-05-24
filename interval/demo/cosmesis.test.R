library(interval)
library(Icens)
data(bcos)
L<-bcos$left
R<-bcos$right
trt<-bcos$treatment

fit1<-icfit(Surv(L,R,type="interval2")~trt)
legend.list<-plot(fit1)
#do.call("legend",legend.list)
#plot(fit1)

summary(fit1)
iR<-EMICM(bcos[bcos$treatment=="Rad",c("left","right")])
iRC<-EMICM(bcos[bcos$treatment=="RadChem",c("left","right")])
class(iR)<-"icfit"
class(iRC)<-"icfit"
summary(iR)
summary(iRC)


fit<-icfit(L,R)

all.scores<-data.frame(L=L,R=R,logrank1=wlr_trafo(L,R,scores="logrank1",icFIT=fit),
   logrank2=wlr_trafo(L,R,scores="logrank2",icFIT=fit),
   wmw=wlr_trafo(L,R,scores="wmw",icFIT=fit))
all.scores

## show that both logrank scores are very close to each other 
plot(all.scores[,3],all.scores[,4],xlab="logrank1",ylab="logrank2",xlim=c(-3,2),ylim=c(-3,2))
lines(c(-10,10),c(-10,10),col=grey(.5))
## show that the wmw are different from the logrank1 scores
plot(all.scores[,3],all.scores[,5],xlab="logrank1",ylab="wmw",xlim=c(-3,2),ylim=c(-3,2))
lines(c(-10,10),c(-10,10),col=grey(.5))

# Network algorithm does not really work on this large of a data set
# use exact.mc to estimate exact p-values by Monte Carlo
#t0<-proc.time()
#out.exact<-ictest(L,R,trt,icFIT=fit,exact=TRUE,mcontrol=mControl(nmc=10^3-1))
#out.exact
#t1<-proc.time()
#t1-t0

## check against results in Fay, 1996, Biometrics 52:811-822
## Note icFIT=fit is not needed just increases speed
ictest(L,R,trt,icFIT=fit,method="pclt",scores="wmw")
ictest(L,R,trt,icFIT=fit,method="pclt",scores="logrank2")
ictest(L,R,trt,icFIT=fit,method="scoretest",alternative="less",scores="wmw")
ictest(L,R,trt,icFIT=fit,method="scoretest",scores="logrank2")



ictest(L,R,trt,icFIT=fit,method="wsr.HLY",alternative="less",scores="wmw")
ictest(L,R,trt,icFIT=fit,method="wsr.HLY",scores="logrank1")
try(ictest(L,R,trt,icFIT=fit,method="wsr.HLY",scores="logrank2"))
ictest(L,R,trt,icFIT=fit,method="wsr.pclt",alternative="less",scores="wmw")
ictest(L,R,trt,icFIT=fit,method="wsr.pclt",scores="logrank1")
try(ictest(L,R,trt,icFIT=fit,method="wsr.pclt",scores="logrank2"))


ictest(L,R,trt,icFIT=fit,method="wsr.mc",scores="wmw",alternative="less",mcontrol=mControl(nwsr=99,np=99))
ictest(L,R,trt,icFIT=fit,method="wsr.mc",scores="logrank1",mcontrol=mControl(nwsr=99,np=99))
ictest(L,R,trt,icFIT=fit,method="wsr.mc",scores="logrank2",mcontrol=mControl(nwsr=99,np=99))

ictest(L,R,trt,icFIT=fit,method="wsr.mc",scores="logrank1")
try(ictest(L,R,trt,icFIT=fit,method="wsr.mc",scores="logrank2"))



## check exact methods, a subset of the data
I<-1:10
tmptrt<-c(rep("Rad",5),rep("RadChem",5))
tmptrt
fitI<-icfit(L[I],R[I])
ictest(L[I],R[I],tmptrt,icFIT=fitI,method="exact.network",alternative="less",scores="wmw")
ictest(L[I],R[I],tmptrt,icFIT=fitI,method="exact.ce",alternative="less",scores="wmw")
ictest(L[I],R[I],tmptrt,icFIT=fitI,method="exact.mc",alternative="less",scores="wmw")
dim(bcos)

## check subsetting
set.seed(1)
I<-rbinom(94,1,.5)==1
icout<-ictest(Surv(left,right,type="interval2")~trt,data=bcos,subset=I)
ictest(bcos$left[I],bcos$right[I],trt[I],initfit=icout$fit)

