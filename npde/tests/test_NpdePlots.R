library(npde)
cat("Running tests on plots\n")

########### Running npde

data(theopp)
data(simtheopp)

xbase1<-autonpde(namobs=theopp,namsim=simtheopp,boolsave=FALSE)

dat3<-theopp
vec<-rep(1,dim(dat3)[1])
loq<-2
vec<-dat3$Wt
for(i in 2:length(vec)) if(is.na(vec[i])) vec[i]<-vec[(i-1)]
dat3$Wt<-vec
vec[!is.na(dat3$Conc) & dat3$Conc>=loq]<-0
dat3$Conc[!is.na(dat3$Conc) & vec==1]<-loq
sex<-ifelse(dat3$Wt>70,"M","F")
dat3<-cbind(dat3,Sex=sex,MDV=rep(0,dim(dat3)[1]),cens=vec)

dat4<-dat3
dat4$ipred<-dat4$Conc*(1+rnorm(dim(dat4)[1],0,0.1))

# Removing censored data
obj0<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="cens1",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="omit")

##### Imputing pd for censored data
## default method (cdf)
obj1<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="cens2",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="cdf")

## ipred
obj2<-autonpde(namobs=dat4,namsim=simtheopp,iid=1,ix=3,iy=4,iipred=9, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="cens2",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="ipred")

## ypred
obj3<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="cens2",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="ypred")

## loq
obj4<-autonpde(namobs=dat3,namsim=simtheopp,iid=1,ix=3,iy=4, imdv=7,icens=8,icov=c(5,6),boolsave=TRUE,type.graph="eps",namsav="cens2",calc.pd=TRUE,calc.npde=TRUE,verbose=FALSE,decorr.method="cholesky",cens.method="loq")

########### Plots

cat("Plots...\n")
# No censoring
plot(xbase1)
plot(xbase1,plot.type="vpc",main="VPC of uncensored data",col="DarkGreen")

# Observed data
par(mfrow=c(2,2))
plot(obj4,plot.type="data",impute.loq=F,new=F, main="Censoring method: loq")
plot(obj2,plot.type="data",impute.loq=F,new=F, main="Censoring method: ipred")
plot(obj3,plot.type="data",impute.loq=F,new=F, main="Censoring method: ypred")
plot(obj1,plot.type="data",impute.loq=F,new=F, main="Censoring method: cdf")
par(mfrow=c(1,1))
title(sub="Raw data, BQL plotted at LOQ",col.sub="Red")

# Complete data
par(mfrow=c(2,2))
plot(obj4,plot.type="data",impute.loq=T,new=F, main="Censoring method: loq")
plot(obj2,plot.type="data",impute.loq=T,new=F, main="Censoring method: ipred")
plot(obj3,plot.type="data",impute.loq=T,new=F, main="Censoring method: ypred")
plot(obj1,plot.type="data",impute.loq=T,new=F, main="Censoring method: cdf")
par(mfrow=c(1,1))
title(sub="Raw data, BQL imputed",col.sub="Red")

# pd, npd
plot(obj1,which="pd")
par(mfrow=c(1,1))
title(sub="Original data (uncensored), pd",col.sub="Red")

plot(obj1,which="pd")
par(mfrow=c(1,1))
title(sub="Method cdf, pd",col.sub="Red")

plot(obj1,which="npd")
par(mfrow=c(1,1))
title(sub="Method cdf, npd",col.sub="Red")

par(mfrow=c(1,2))
plot(obj1,which="npd",plot.type="ecdf",new=F, main="Censoring method: cdf")
plot(obj1,which="pd",plot.type="ecdf",new=F)
par(mfrow=c(1,1))
title(sub="Method cdf, ecdf for npd & pd",col.sub="Red")

par(mfrow=c(2,2))
plot(obj3,which="pd",plot.type="qqplot",new=F, main="Censoring method: omit")
plot(obj4,which="pd",plot.type="qqplot",new=F, main="Censoring method: loq")
plot(obj2,which="pd",plot.type="qqplot",new=F, main="Censoring method: ipred")
plot(obj1,which="pd",plot.type="qqplot",new=F, main="Censoring method: cdf")
par(mfrow=c(1,1))
title(sub="QQplots",col.sub="Red")

# npde
plot(obj1)
par(mfrow=c(1,1))
title(sub="Default plots, original data (no cens)",col.sub="Red")

plot(obj1)
par(mfrow=c(1,1))
title(sub="Default plots, cens=cdf",col.sub="Red")

plot(obj2)
par(mfrow=c(1,1))
title(sub="Default plots, cens=ipred",col.sub="Red")

par(mfrow=c(1,1))
plot(obj1,plot.type="x.scatter",new=F, main="Censoring method: cdf")
title(sub="Scatter-plot, cens=cdf",col.sub="Red")

# VPC
plot(obj0,plot.type="vpc")
title(sub="VPC, original data (no cens)",col.sub="Red")

plot(obj3,plot.type="vpc")
title(sub="VPC, cens=omit, BQL data omitted (default)",col.sub="Red")

plot(obj3,plot.type="vpc",impute.loq=TRUE)
title(sub="VPC, cens=omit, BQL data plotted at LOQ",col.sub="Red")

plot(obj1,plot.type="vpc",impute.loq=FALSE)
title(sub="VPC, cens=cdf, BQL data plotted at LOQ",col.sub="Red")

plot(obj1,plot.type="vpc")
title(sub="VPC, cens=cdf, BQL data imputed (default)",col.sub="Red")

plot(obj2,plot.type="vpc")
title(sub="VPC, cens=ipred, BQL data imputed (default)",col.sub="Red")

plot(obj1,plot.type="vpc",ylog=TRUE)
title(sub="VPC, cens=cdf, BQL data imputed (default), log-scale",col.sub="Red")

# P(LOQ)
plot(obj1,plot.type="loq")
title(sub="P(y<LOQ), cens=cdf",col.sub="Red")

# Split by covariates
par(mfrow=c(2,3))
plot(obj1,which="npd",plot.type="ecdf",covsplit=T,new=F)
plot(1:10,1:10,xlab="",ylab="",axes=F,type="n")
text(5,9,"Covariate plots",col="Red",cex=2)
text(5,7,"npd emp cdf",col="Red",cex=2)
text(5,5,"cens=cdf",col="Red",cex=2)

par(mfrow=c(2,3))
plot(obj1,which="npde",plot.type="x.scatter",covsplit=T,new=F)
plot(1:10,1:10,xlab="",ylab="",axes=F,type="n")
text(5,9,"Covariate plots",col="Red",cex=2)
text(5,7,"npde vs X",col="Red",cex=2)
text(5,5,"cens=cdf",col="Red",cex=2)

cat("End of tests on plots...\n")

####################################################################################
####				Test for interactive menu				####
####################################################################################

# Menu
# source("../R/func_menus.R")
# xinput<-pdemenu()

####################################################################################
