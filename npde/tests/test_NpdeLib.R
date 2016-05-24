library(npde)
cat("Running tests for the main computational functions\n")

data(theopp)
data(simtheopp)

xbase1<-autonpde(namobs=theopp,namsim=simtheopp,boolsave=FALSE)

# Cholesky decomposition
xbase1<-autonpde(namobs=theopp,namsim=simtheopp, iid=1,ix=3,iy=4,imdv=0,icens=0,icov=0,boolsave=TRUE,type.graph="eps",
namsav="base1",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky", cens.method="omit",units=list(x="hr",y="mg/L"))

# Decorrelation through the inverse method
xbase2<-autonpde(namobs=theopp,namsim=simtheopp, iid=1,ix=3,iy=4,imdv=0,icens=0,icov=0,boolsave=TRUE,type.graph="eps",
namsav="base2",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="inverse", cens.method="omit",units=list(x="hr",y="mg/L"))

idx<-xbase1@data@not.miss
ytest1<-max(abs(xbase1@results@res$pd[idx]-xbase2@results@res$pd[idx]))

if(zapsmall(ytest1)==0) cat("Cholesky/inverse method: same pd, test successful\n") else cat("Problem with xbase1/xbase2\n")

####### With censoring

dat3<-theopp
vec<-rep(1,dim(dat3)[1])
loq<-2
vec[!is.na(dat3$Conc) & dat3$Conc>=loq]<-0
dat3$Conc[!is.na(dat3$Conc) & vec==1]<-loq
sex<-ifelse(dat3$Wt>70,"M","F")
dat3<-cbind(dat3,Sex=sex,MDV=rep(0,dim(dat3)[1]),cens=vec)

dat4<-dat3
dat4$ipred<-dat4$Conc*(1+rnorm(dim(dat4)[1],0,0.1))

# Removing censored data
xcens1<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="cens1",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="omit")

##### Imputing pd for censored data
## default method (cdf)
xcens2<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="cens2",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="cdf")

## method ipred (automatic detection of ipred)
xipred<-autonpde(namobs=dat4,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps", namsav="ipred1",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="ipred")

## method loq
xloq<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="loq1",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="loq")

## method ypred
xypred<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="ypred1",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="ypred")

# Automatic column recognition
# xdetect<-autonpde(namobs=dat4,namsim=simtheopp,icov=c(5,6), namsav="idetect"))

# Without all defaults
xcens1.bis<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),namsav="cens1.bis",cens.method="omit")

xcens2.bis<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),namsav="cens2.bis")

idx<-(dat3$MDV==0 & dat3$cens==0)
ytest1<-max(abs(xcens1@results@res$pd[idx]-xcens1.bis@results@res$pd[idx]))

ytest2<-max(abs(xcens2@results@res$pd[idx]-xcens2.bis@results@res$pd[idx]))

ytest3<-max(abs(xcens2@results@res$pd[idx]-xipred@results@res$pd[idx]))

if(zapsmall(ytest1+ytest2+ytest3)==0) cat("Test with censoring successful \n") else cat("Different results, check defaults\n")

cat("Testing gof.test function\n")

gof.test(xbase1)

gof.test(xcens2)

ygof<-gof.test(xcens2,which="pd",covsplit=TRUE)
print(ygof)

gof.test(xcens2,which="pd",sample=TRUE)

cat("End of tests for the main computational functions\n")

