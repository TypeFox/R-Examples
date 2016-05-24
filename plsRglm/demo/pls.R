set.seed(12345)

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
plsR(yCornell,XCornell,10)$InfCrit
plsR(yCornell,XCornell,10,typeVC="standard")$InfCrit
plsR(yCornell,XCornell,6)$AIC 
plsR(yCornell,XCornell,6)$AIC.std    
plsRglm(yCornell,XCornell,3)$uscores
plsRglm(yCornell,XCornell,3)$pp
plsRglm(yCornell,XCornell,3)$Coeffs
plsRglm(yCornell,XCornell,10)$InfCrit
plsRglm(yCornell,XCornell,10,modele="pls-glm-gaussian")$InfCrit


bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb,MClassed=TRUE)
bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,keepcoeffs=TRUE)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb,MClassed=TRUE)

bbb <- PLS_glm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls")
kfolds2CVinfos_glm(bbb)


PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=12,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=12,keepfolds=FALSE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=6,NK=2,random=FALSE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=6,NK=2,random=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE,keepfolds=FALSE)

bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=6,NK=1)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb2)
PLS_lm(yCornell,XCornell,6,typeVC="standard")$InfCrit


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

modpls <- PLS_lm(yaze_compl,Xaze_compl,10,MClassed=TRUE)

modpls$AIC
modpls$AIC.std
modpls$MissClassed
modpls$Probs
modpls$Probs.trc
modpls$Probs-modpls$Probs.trc

modpls$InfCrit
PLS_lm(yaze_compl,Xaze_compl,10)$InfCrit

PLS_lm(yaze_compl,Xaze_compl,10,typeVC="standard")$InfCrit
PLS_lm(yaze_compl,Xaze_compl,10,typeVC="standard",MClassed=TRUE)$InfCrit
rm(list=c("Xaze_compl","yaze_compl","modpls"))


dimX <- 24
Astar <- 2
simul_data_UniYX(dimX,Astar)
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ydataAstar2 <- dataAstar2[,1]
XdataAstar2 <- dataAstar2[,2:(dimX+1)]
ysimbin1 <- dicho(ydataAstar2)
Xsimbin1 <- dicho(XdataAstar2)
PLS_lm(ysimbin1,Xsimbin1,10,MClassed=TRUE)$Probs
PLS_lm(ysimbin1,Xsimbin1,10,MClassed=TRUE)$Probs.trc
PLS_lm(ysimbin1,Xsimbin1,10,MClassed=TRUE)$MissClassed
PLS_lm(ysimbin1,Xsimbin1,10,typeVC="standard",MClassed=TRUE)$InfCrit
PLS_lm(ysimbin1,XdataAstar2,10,typeVC="standard",MClassed=TRUE)$InfCrit
PLS_lm(ydataAstar2,XdataAstar2,10,typeVC="standard")$InfCrit
rm(list=c("dimX","Astar","dataAstar2","ysimbin1","Xsimbin1","ydataAstar2","XdataAstar2"))


dimX <- 6
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ydataAstar4 <- dataAstar4[,1]
XdataAstar4 <- dataAstar4[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar4,XdataAstar4,10,typeVC="standard")
modpls$computed_nt
modpls$InfCrit
str(modpls)
rm(list=c("dimX","Astar","dataAstar4","modpls","ydataAstar4","XdataAstar4"))


dimX <- 24
Astar <- 2
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ydataAstar2 <- dataAstar2[,1]
XdataAstar2 <- dataAstar2[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar2,XdataAstar2,10,typeVC="standard")
modpls$computed_nt
modpls$InfCrit
rm(list=c("dimX","Astar","dataAstar2","modpls","ydataAstar2","XdataAstar2"))


data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

set.seed(1385)
Cornell.tilt.boot <- tilt.bootpls(plsR(yCornell,XCornell,3), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)
Cornell.tilt.boot
str(Cornell.tilt.boot)

boxplots.bootpls(Cornell.tilt.boot,indices=2:7)

rm(Cornell.tilt.boot)


# Comparing the results with the plspm package and SIMCA results in Tenenhaus's book
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_lm(yCornell,XCornell,3)$uscores
PLS_lm(yCornell,XCornell,3)$pp
PLS_lm(yCornell,XCornell,3)$Coeffs

PLS_lm(yCornell,XCornell,4,typeVC="standard")$press.ind
PLS_lm(yCornell,XCornell,4,typeVC="standard")$press.tot
PLS_lm(yCornell,XCornell,4,typeVC="standard")$InfCrit

PLS_lm_wvc(dataY=yCornell,dataX=XCornell,nt=3,dataPredictY=XCornell[1,])
PLS_lm_wvc(dataY=yCornell[-c(1,2)],dataX=XCornell[-c(1,2),],nt=3,dataPredictY=XCornell[c(1,2),])
PLS_lm_wvc(dataY=yCornell[-c(1,2)],dataX=XCornell[-c(1,2),],nt=3,dataPredictY=XCornell[c(1,2),],keepcoeffs=TRUE)


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
PLS_lm(log(ypine),Xpine,4)$Std.Coeffs
PLS_lm(log(ypine),Xpine,4)$Coeffs
PLS_lm(log(ypine),Xpine,1)$Std.Coeffs
PLS_lm(log(ypine),Xpine,1)$Coeffs
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$InfCrit


data(pine_full)
Xpine_full<-pine_full[,1:10]
ypine_full<-pine_full[,11]
PLS_lm(log(ypine_full),Xpine_full,1)$Std.Coeffs
PLS_lm(log(ypine_full),Xpine_full,1)$Coeffs
cor(cbind(Xpine,log(ypine)))


XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
PLS_lm(log(ypine),XpineNAX21,4)$Std.Coeffs

PLS_lm(log(ypine),XpineNAX21,4)$YChapeau[1,]
PLS_lm(log(ypine),Xpine,4)$YChapeau[1,]

PLS_lm(log(ypine),XpineNAX21,4)$CoeffC

PLS_lm(log(ypine),XpineNAX21,2,dataPredictY=XpineNAX21[1,])$ValsPredictY

PLS_lm(log(ypine),Xpine,10,typeVC="none")$InfCrit
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$InfCrit
PLS_lm(log(ypine),Xpine,10,typeVC="adaptative")$InfCrit
PLS_lm(log(ypine),Xpine,10,typeVC="missingdata")$InfCrit
PLS_lm(log(ypine),XpineNAX21,10,typeVC="none")$InfCrit
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$InfCrit
PLS_lm(log(ypine),XpineNAX21,10,typeVC="adaptative")$InfCrit
PLS_lm(log(ypine),XpineNAX21,10,typeVC="missingdata")$InfCrit

PLS_lm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeau
PLS_lm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeauNA

# The results from plspm were uncorrect
if ("plspm" %in%  installed.packages()){
library(plspm)
plsreg1(x=XCornell,y=as.vector(yCornell),nc=3)$coeffs
plsreg1(x=XCornell,y=as.vector(yCornell),nc=4,cv=TRUE)
plsreg1(x=XCornell,y=as.vector(yCornell),nc=4,cv=TRUE)$Q2
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=4)$std.coef
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=4)$coeffs
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=10)$Q2
plsreg1(x=Xpine,y=log(as.vector(ypine)),nc=4,cv=TRUE)$Q2
plsreg1(x=Xpine_full,y=log(as.vector(ypine_full)),nc=4)$coeffs
plsreg1(x=XpineNAX21,y=as.vector(log(ypine)),nc=4,cv=TRUE)
}


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=Xpine,nt=10,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=Xpine,nt=10,K=6,NK=1)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb2)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$InfCrit

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbbNA <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=XpineNAX21,nt=10,K=12,NK=1)
bbbNA2 <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=XpineNAX21,nt=10,K=6,NK=1)
kfolds2CVinfos_lm(bbbNA)
kfolds2CVinfos_lm(bbbNA2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$InfCrit


data(XpineNAX21)
PLS_lm_wvc(dataY=log(ypine)[-1],dataX=XpineNAX21[-1,],nt=3)
PLS_lm_wvc(dataY=log(ypine)[-1],dataX=XpineNAX21[-1,],nt=3,dataPredictY=XpineNAX21[1,])
PLS_lm_wvc(dataY=log(ypine)[-2],dataX=XpineNAX21[-2,],nt=3,dataPredictY=XpineNAX21[2,])
PLS_lm_wvc(dataY=log(ypine),dataX=XpineNAX21,nt=3)
rm(list=c("Xpine","XpineNAX21","ypine","bbb","bbb2","bbbNA","bbbNA2"))


dimX <- 24
Astar <- 3
dataAstar3 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar3 <- dataAstar3[,1]
XdataAstar3 <- dataAstar3[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar3,XdataAstar3,10,typeVC="standard")
modpls$computed_nt
modpls$InfCrit
rm(list=c("dimX","Astar","dataAstar3","modpls","ydataAstar3","XdataAstar3"))


dimX <- 24
Astar <- 4
dataAstar4 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar4 <- dataAstar4[,1]
XdataAstar4 <- dataAstar4[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar4,XdataAstar4,10,typeVC="standard")
modpls$computed_nt
modpls$InfCrit
rm(list=c("dimX","Astar","dataAstar4","modpls","ydataAstar4","XdataAstar4"))


dimX <- 24
Astar <- 5
dataAstar5 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar5 <- dataAstar5[,1]
XdataAstar5 <- dataAstar5[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar5,XdataAstar5,10,typeVC="standard")
modpls$computed_nt
modpls$InfCrit
rm(list=c("dimX","Astar","dataAstar5","modpls","ydataAstar5","XdataAstar5"))


dimX <- 24
Astar <- 6
dataAstar6 <- t(replicate(200,simul_data_UniYX(dimX,Astar)))
ydataAstar6 <- dataAstar6[,1]
XdataAstar6 <- dataAstar6[,2:(dimX+1)]
modpls <- PLS_lm(ydataAstar6,XdataAstar6,10,typeVC="standard")
modpls$computed_nt
modpls$InfCrit
rm(list=c("dimX","Astar","dataAstar6","modpls","ydataAstar6","XdataAstar6"))




# Lazraq-Cleroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)
boot::boot.array(Cornell.boot, indices=TRUE)

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.boot$t[,-1])~factor(rep(1:7,rep(250,7))), main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.boot$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.boot,indices=2:8)
# Confidence intervals plotting
confints.bootpls(Cornell.boot,indices=2:8)
plots.confints.bootpls(confints.bootpls(Cornell.boot,indices=2:8))


data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cleroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)

(temp.ci <- confints.bootpls(Cornell.boot,2:8))
plots.confints.bootpls(temp.ci)
(temp.ci <- confints.bootpls(Cornell.boot,2:8,typeBCa=FALSE))
plots.confints.bootpls(temp.ci)
(temp.ci <- confints.bootpls(Cornell.boot,c(2,4,6)))
plots.confints.bootpls(temp.ci)


plot(Cornell.boot,index=2)
boot::jack.after.boot(Cornell.boot, index=2, useJ=TRUE, nt=3)
plot(Cornell.boot,index=2,jack=TRUE)

car::dataEllipse(Cornell.boot$t[,2], Cornell.boot$t[,3], cex=.3, levels=c(.5, .95, .99), robust=T)


rm(Cornell.boot)


# PLS balanced bootstrap

set.seed(225)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="balanced", stype="i", R=250)
boot::boot.array(Cornell.boot, indices=TRUE)



# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.boot$t[,-1])~factor(rep(1:7,rep(250,7))), main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.boot$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.boot,indices=2:8)
# Confidence intervals plotting
confints.bootpls(Cornell.boot,indices=2:8)
plots.confints.bootpls(confints.bootpls(Cornell.boot,indices=2:8))


library(boot)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=7)
boot::boot.ci(Cornell.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=8)

plot(Cornell.boot,index=2)
boot::jack.after.boot(Cornell.boot, index=2, useJ=TRUE, nt=3)
plot(Cornell.boot,index=2,jack=TRUE)

rm(Cornell.boot)

# PLS permutation bootstrap

set.seed(500)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="permutation", stype="i", R=1000)
boot::boot.array(Cornell.boot, indices=TRUE)


# Graph of bootstrap distributions
boxplot(as.vector(Cornell.boot$t[,-1])~factor(rep(1:7,rep(1000,7))),main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.boot$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.boot,indices=2:8)




library(boot)
plot(Cornell.boot,index=2)

qqnorm(Cornell.boot$t[,2],ylim=c(-1,1))
abline(h=Cornell.boot$t0[2],lty=2)
(sum(abs(Cornell.boot$t[,2])>=abs(Cornell.boot$t0[2]))+1)/(length(Cornell.boot$t[,2])+1)

rm(Cornell.boot)


























data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
plsRglm(log(ypine),Xpine,1)$Std.Coeffs
plsRglm(log(ypine),Xpine,1)$Coeffs
plsRglm(log(ypine),Xpine,4)$Std.Coeffs
plsRglm(log(ypine),Xpine,4)$Coeffs
plsRglm(log(ypine),Xpine,4)$PredictY[1,]
plsRglm(log(ypine),Xpine,4,dataPredictY=Xpine[1,])$PredictY[1,]

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
str(plsRglm(log(ypine),XpineNAX21,2))
plsRglm(log(ypine),XpineNAX21,4)$Std.Coeffs
plsRglm(log(ypine),XpineNAX21,4)$YChapeau[1,]
plsRglm(log(ypine),Xpine,4)$YChapeau[1,]
plsRglm(log(ypine),XpineNAX21,4)$CoeffC
plsRglm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeau
plsRglm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeauNA




