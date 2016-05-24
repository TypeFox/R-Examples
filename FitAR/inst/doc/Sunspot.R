#Compare spectral density function for monthly sunspots numbers
library(FitAR) #needs Version 1.78 or higher
#First model order selection
pARBIC<-SelectModel(sunspots, lag.max=300, ARModel="AR", Criterion="BIC", Best=1)
pUBIC<-SelectModel(sunspots, lag.max=300, ARModel="ARz", Criterion="UBIC", Best=1)
pBIC<-SelectModel(sunspots, lag.max=300, ARModel="ARz", Criterion="BIC", Best=1)
pQBIC02<-SelectModel(sunspots, lag.max=300, ARModel="ARz", Criterion="QBIC", Best=1, Q=0.02)
#
#Fit models
ansARBIC<-FitAR(sunspots, p=pARBIC, ARModel="ARp")
ansBIC<-FitAR(sunspots, p=pBIC, ARModel="ARz")
ansUBIC<-FitAR(sunspots, p=pUBIC, ARModel="ARz")
ansQBIC02<-FitAR(sunspots, p=pQBIC02, ARModel="ARz")
#
#multi-panel plot, standard graphics
par(mfrow=c(2,2))
sARBIC<-sdfplot(ansARBIC)
title(main="AR, BIC")
sBIC<-sdfplot(ansBIC)
title(main="Subset, BIC")
sUBIC<-sdfplot(ansUBIC)
title(main="Subset, UBIC")
sQBIC<-sdfplot(ansQBIC02)
title(main="Subset, BICq(q=0.02)")
par(mfrow=c(1,1))
#
#improved lattice plot
n<-nrow(sARBIC)
which<-ordered(rep(c("AR/BIC","Subset/BIC","Subset/UBIC","Subset/QBIC"),rep(n,4)))
SDF<-rbind(sARBIC,sBIC,sUBIC,sQBIC)
SDF[,1]<-SDF[,1]/pi
S.df <- data.frame(f=SDF[,1], sdf=SDF[,2], which=which)
lattice.options(par.strip.text=list(col="white"))
trellis.device(color=FALSE) #no color for print publication
xyplot(sdf~f|which, type="l", xlab="frequency", ylab="log spectral density", data=S.df)
#
#Include AIC: full AR and subset AR
#pure ar
pARAIC<-SelectModel(sunspots, lag.max=300, ARModel="AR", Criterion="AIC", Best=1)
ansARAIC<-FitAR(sunspots, p=pARAIC, ARModel="ARp")
sARAIC<-sdfplot(ansARAIC)
#subset
pAIC<-SelectModel(sunspots, lag.max=300, ARModel="ARz", Criterion="AIC", Best=1)
ansAIC<-FitAR(sunspots, p=pAIC, ARModel="ARz")
sAIC<-sdfplot(ansAIC)
#
NamesModels<-c("AR/AIC","AR/BIC","Subset/AIC", "Subset/BIC","Subset/UBIC","Subset/QBIC")
whichModels<-ordered(x=rep(NamesModels,rep(n,6)),levels=NamesModels)
SDF<-rbind(sARAIC, sARBIC, sAIC, sBIC, sUBIC, sQBIC)
SDF[,1]<-SDF[,1]/pi
S.df <- data.frame(f=SDF[,1], sdf=SDF[,2], which=whichModels)
lattice.options(par.strip.text=list(col="white"))
trellis.device(color=FALSE) #no color for print publication
xyplot(sdf~f|which, type="l", xlab="frequency", ylab="log spectral density", data=S.df)
#
#Compare the deviances
devTable<-matrix(numeric(2*length(NamesModels)), nrow=2)
devTable[1,]<--2*c(ansARAIC$loglikelihood,ansARBIC$loglikelihood,ansAIC$loglikelihood,ansBIC$loglikelihood,ansUBIC$loglikelihood,ansQBIC02$loglikelihood)
devTable[2,]<-c(length(ansARAIC$pvec),length(ansARBIC$pvec),length(ansAIC$pvec),length(ansBIC$pvec),length(ansUBIC$pvec),length(ansQBIC02$pvec))
dimnames(devTable)<-list(c("deviance","NumberParameters"), NamesModels)
devTable


#find q so BICq equivalent to UBIC - use trial-and-error method
#answer: (0.02, 0.21)
#This script takes a few minutes
lbic<-length(SelectModel(sunspots, lag.max=300, ARModel="ARz", Criterion="UBIC", Best=1))
QS<-seq(0.02, 0.21, 0.01)
wq<-logical(length(QS))
sz<-numeric(length(QS))
for (i in 1:length(QS)){
    qS<-length(SelectModel(sunspots, lag.max=300, ARModel="ARz", Criterion="QBIC", Best=1, Q=QS[i]))
    sz[i] <- qS
    wq[i]<-lbic==qS
    }
c(min(QS[wq]),max(QS[wq]))
