set.seed(12345)

data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
PLS_glm(log(ypine),Xpine,1)$Std.Coeffs
PLS_glm(log(ypine),Xpine,1)$Coeffs
PLS_glm(log(ypine),Xpine,4)$Std.Coeffs
PLS_glm(log(ypine),Xpine,4)$Coeffs
PLS_glm(log(ypine),Xpine,4)$PredictY[1,]
PLS_glm(log(ypine),Xpine,4,dataPredictY=Xpine[1,])$PredictY[1,]

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
str(PLS_glm(log(ypine),XpineNAX21,2))
PLS_glm(log(ypine),XpineNAX21,4)$Std.Coeffs
PLS_glm(log(ypine),XpineNAX21,4)$YChapeau[1,]
PLS_glm(log(ypine),Xpine,4)$YChapeau[1,]
PLS_glm(log(ypine),XpineNAX21,4)$CoeffC
PLS_glm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeau
PLS_glm(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeauNA

# compare pls-glm-gaussian with classic plsR
cbind(PLS_glm(log(ypine),Xpine,4,modele="pls")$Std.Coeffs,PLS_glm(log(ypine),Xpine,4,modele="pls-glm-gaussian")$Std.Coeffs)

# without missing data
cbind(log(ypine),PLS_glm(log(ypine),Xpine,4,modele="pls")$YChapeau,PLS_glm(log(ypine),Xpine,4,modele="pls-glm-gaussian")$YChapeau)
cbind(log(ypine),PLS_glm(log(ypine),XpineNAX21,4,modele="pls")$YChapeau,PLS_glm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$YChapeau)

# with missing data
cbind((log(ypine)),PLS_glm(log(ypine),XpineNAX21,4,modele="pls")$YChapeau,PLS_glm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$YChapeau)
cbind((log(ypine)),PLS_glm(log(ypine),XpineNAX21,4,modele="pls")$ValsPredictY,PLS_glm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$ValsPredictY)


# compare pls-glm-gaussian with log link with classic plsR on the log
cbind(PLS_glm(log(ypine),Xpine,4,modele="pls")$Std.Coeffs,PLS_glm(log(ypine),Xpine,4,modele="pls-glm-gaussian")$Std.Coeffs,PLS_glm(log(ypine),Xpine,4,modele="pls-glm-family",family=gaussian(link="identity"))$Std.Coeffs,PLS_glm(ypine,Xpine,4,modele="pls-glm-family",family=gaussian(link=log))$Std.Coeffs)

# without missing data
cbind(log(ypine),PLS_glm(log(ypine),Xpine,4,modele="pls")$YChapeau,PLS_glm(log(ypine),Xpine,4,modele="pls-glm-gaussian")$YChapeau,PLS_glm(log(ypine),Xpine,4,modele="pls-glm-family",family=gaussian(link="identity"))$YChapeau,log(PLS_glm(ypine,XpineNAX21,4,modele="pls-glm-family",family=gaussian(link=log))$YChapeau))

# with missing data
cbind((log(ypine)),PLS_glm(log(ypine),XpineNAX21,4,modele="pls")$YChapeau,PLS_glm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$YChapeau,PLS_glm(log(ypine),XpineNAX21,4,modele="pls-glm-family",family=gaussian(link="identity"))$YChapeau,log(PLS_glm(ypine,XpineNAX21,4,modele="pls-glm-family",family=gaussian(link=log))$YChapeau))
cbind((log(ypine)),PLS_glm(log(ypine),XpineNAX21,4,modele="pls")$ValsPredictY,PLS_glm(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$ValsPredictY,PLS_glm(log(ypine),XpineNAX21,4,modele="pls-glm-family",family=gaussian(link="identity"))$ValsPredictY,log(PLS_glm(ypine,XpineNAX21,4,modele="pls-glm-family",family=gaussian(link=log))$ValsPredictY))


#other links
data.frame(pls=PLS_glm(log(ypine),Xpine,4,modele="pls")$Std.Coeffs,identity=PLS_glm(log(ypine),Xpine,4,modele="pls-glm-family",family=gaussian(link="identity"))$Std.Coeffs,log=PLS_glm(ypine,Xpine,4,modele="pls-glm-family",family=gaussian(link=log))$Std.Coeffs,inverse=PLS_glm(ypine,Xpine,4,modele="pls-glm-family",family=gaussian(link=inverse))$Std.Coeffs)


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
modpls$valpvalstep
modpls$Coeffsmodel_vals
modpls2 <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-family",family=binomial(link=logit),MClassed=TRUE,pvals.expli=TRUE)
modpls3 <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-family",family=binomial(link=probit),MClassed=TRUE,pvals.expli=TRUE)
modpls4 <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-family",family=binomial(link=cauchit),MClassed=TRUE,pvals.expli=TRUE)
#fails modpls5 <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-family",family=binomial(link=log),MClassed=TRUE,pvals.expli=TRUE)
modpls6 <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-family",family=binomial(link=cloglog),MClassed=TRUE,pvals.expli=TRUE)
data.frame(logit=modpls2$Std.Coeffs,probit=modpls3$Std.Coeffs,cauchit=modpls4$Std.Coeffs,log=NA,cloglog=modpls6$Std.Coeffs)




# Test of other families and links on the same datasets.
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
modpls <- PLS_glm(ypine,Xpine,10,modele="pls")
modpls$computed_nt
modpls$InfCrit
modpls2 <- PLS_glm(ypine,Xpine,10,modele="pls-glm-Gamma")
modpls2$InfCrit
modpls2a <- PLS_glm(ypine,Xpine,10,modele="pls-glm-family",family=Gamma(link=inverse))
modpls2a$InfCrit
modpls2b <- PLS_glm(ypine+1,Xpine,10,modele="pls-glm-family",family=Gamma(link=identity))
modpls2b$InfCrit
modpls2c <- PLS_glm(ypine,Xpine,10,modele="pls-glm-family",family=Gamma(link=log))
modpls2c$InfCrit


modpls3 <- PLS_glm(ypine,Xpine,10,modele="pls-glm-gaussian")
modpls3$InfCrit
modpls3a <- PLS_glm(ypine,Xpine,10,modele="pls-glm-family",family=gaussian(link=identity))
modpls3a$InfCrit
modpls3b <- PLS_glm(ypine,Xpine,10,modele="pls-glm-family",family=gaussian(link=log))
modpls3b$InfCrit
modpls3c <- PLS_glm(ypine,Xpine,10,modele="pls-glm-family",family=gaussian(link=inverse))
modpls3c$InfCrit


modpls4 <- PLS_glm(round(ypine),Xpine,10,modele="pls-glm-poisson")
modpls4$InfCrit
modpls4a <- PLS_glm(round(ypine),Xpine,10,modele="pls-glm-family",family=poisson(link=log))
modpls4a$InfCrit
modpls4b <- PLS_glm(round(ypine)+1,Xpine,10,modele="pls-glm-family",family=poisson(link=identity))
modpls4b$InfCrit
#Fails
#modpls4c <- PLS_glm(round(ypine),Xpine,10,modele="pls-glm-family",family=poisson(link=sqrt))
#modpls4c$InfCrit


data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- PLS_glm(yCornell,XCornell,10,modele="pls-glm-inverse.gaussian")
modpls$InfCrit
modplsa <- PLS_glm(yCornell,XCornell,10,modele="pls-glm-family",family=inverse.gaussian(link=1/mu^2))
modplsa$InfCrit
modplsb <- PLS_glm(yCornell,XCornell,10,modele="pls-glm-family",family=inverse.gaussian(link=inverse))
modplsb$InfCrit
modplsc <- PLS_glm(yCornell,XCornell,10,modele="pls-glm-family",family=inverse.gaussian(link=identity))
modplsc$InfCrit
modplsd <- PLS_glm(yCornell,XCornell,10,modele="pls-glm-family",family=inverse.gaussian(link=log))
modplsd$InfCrit






data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_glm_wvc(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",dataPredictY=XCornell[1,])
PLS_glm_wvc(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian(),dataPredictY=XCornell[1,])
PLS_glm_wvc(dataY=yCornell[-1],dataX=XCornell[-1,],nt=3,modele="pls-glm-gaussian",dataPredictY=XCornell[1,])
PLS_glm_wvc(dataY=yCornell[-1],dataX=XCornell[-1,],nt=3,modele="pls-glm-family",family=gaussian(),dataPredictY=XCornell[1,])


## With an incomplete dataset (X[1,2] is NA)
data(pine)
ypine <- pine[,11]
data(XpineNAX21)
PLS_glm_wvc(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian")



data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
PLS_glm_wvc(ypine,Xpine,10,modele="pls")
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-Gamma")
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-family",family=Gamma())
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-gaussian")
PLS_glm_wvc(ypine,Xpine,10,modele="pls-glm-family",family=gaussian(log))
PLS_glm_wvc(round(ypine),Xpine,10,modele="pls-glm-poisson")
PLS_glm_wvc(round(ypine),Xpine,10,modele="pls-glm-family",family=poisson(log))


data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_glm_wvc(yCornell,XCornell,10,modele="pls-glm-inverse.gaussian")
PLS_glm_wvc(yCornell,XCornell,10,modele="pls-glm-family",family=inverse.gaussian())





XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,NK=1,modele="pls")
kfolds2CVinfos_glm(bbb)

PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=12)
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=TRUE,keepfolds=TRUE)$results_kfolds

#Different ways of model specifications
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian,K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds


bbb2 <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,modele="pls-glm-gaussian",keepcoeffs=TRUE)
bbb2 <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),K=6,keepcoeffs=TRUE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(yCornell),XCornell,10,typeVC="standard")$InfCrit


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-family",family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-gaussian",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])

kfolds2Chisqind(bbb)
kfolds2Chisq(bbb)
kfolds2CVinfos_glm(bbb)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$InfCrit

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-family",family=gaussian(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-gaussian",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$InfCrit


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=10,modele="pls",keepcoeffs=TRUE)

#For Jackknife computations
kfolds2coeff(bbb)
bbb2 <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=10,modele="pls-glm-family",family=binomial(probit),keepcoeffs=TRUE)
bbb2 <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=10,modele="pls-glm-logistic",keepcoeffs=TRUE)
kfolds2CVinfos_glm(bbb,MClassed=TRUE)
kfolds2CVinfos_glm(bbb2,MClassed=TRUE)
kfolds2coeff(bbb2)

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)



data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=Xpine,nt=10,modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=Xpine,nt=10,modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])

kfolds2Chisqind(bbb)
kfolds2Chisq(bbb)
kfolds2CVinfos_glm(bbb)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$InfCrit

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb2 <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=XpineNAX21,nt=10,modele="pls-glm-family",family=poisson(log),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb2 <- PLS_glm_kfoldcv(dataY=round(ypine),dataX=XpineNAX21,nt=10,modele="pls-glm-poisson",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$InfCrit



data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-family",family=Gamma,K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb <- PLS_glm_kfoldcv(dataY=ypine,dataX=Xpine,nt=10,modele="pls-glm-Gamma",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])

kfolds2Chisqind(bbb)
kfolds2Chisq(bbb)
kfolds2CVinfos_glm(bbb)
PLS_lm(log(ypine),Xpine,10,typeVC="standard")$InfCrit

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-family",family=Gamma(),K=10,keepcoeffs=TRUE,keepfolds=FALSE)
bbb2 <- PLS_glm_kfoldcv(dataY=ypine,dataX=XpineNAX21,nt=10,modele="pls-glm-Gamma",K=10,keepcoeffs=TRUE,keepfolds=FALSE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(ypine),XpineNAX21,10,typeVC="standard")$InfCrit



data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,NK=1,modele="pls")
kfolds2CVinfos_glm(bbb)

PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=12)
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=inverse.gaussian,K=12)
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,NK=2,random=TRUE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=inverse.gaussian(),K=6,NK=2,random=TRUE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-inverse.gaussian",K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=inverse.gaussian(link = "1/mu^2"),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds

bbb2 <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,modele="pls-glm-inverse.gaussian",keepcoeffs=TRUE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(yCornell),XCornell,10,typeVC="standard")$InfCrit




data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
bbb <- PLS_glm_kfoldcv(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-family",family="binomial")
bbb <- PLS_glm_kfoldcv(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-logistic")
bbb2 <- PLS_glm_kfoldcv(dataY=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-family",family=binomial(),K=10)
bbb2 <- PLS_glm_kfoldcv(dataY=yaze_compl,dataX=Xaze_compl,nt=10,modele="pls-glm-logistic",K=10)
kfolds2Chisq(bbb)
kfolds2Chisqind(bbb)
kfolds2Chisq(bbb2)
kfolds2Chisqind(bbb2)


bbb <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls")
kfolds2CVinfos_glm(bbb,MClassed=TRUE)
bbba <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls-glm-family",family=gaussian())
kfolds2CVinfos_glm(bbba,MClassed=TRUE)
bbb2 <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls-glm-logistic")
kfolds2CVinfos_glm(bbb2,MClassed=TRUE)
bbb2a <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls-glm-family",family=binomial())
kfolds2CVinfos_glm(bbb2a,MClassed=TRUE)


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

dataset <- cbind(y=yaze_compl,Xaze_compl)

library(boot)
# Lazraq-Cleroux PLS bootstrap Classic

#aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), typeboot="plsmodel", sim="ordinary", stype="i", R=250)
# The same
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-family",family=binomial), typeboot="plsmodel", sim="ordinary", stype="i", R=250)

# Bastien CSDA 2005 Bootstrap

#aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), typeboot="fmodel_np", sim="ordinary", stype="i", R=250)
# The same
aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-family",family=binomial), typeboot="fmodel_np", sim="ordinary", stype="i", R=250)





data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

dataset <- cbind(y=yaze_compl,Xaze_compl)

library(boot)
# Lazraq-Cleroux PLS bootstrap Classic

aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic", family=NULL), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), typeboot="plsmodel", sim="ordinary", stype="i", index=1)
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), typeboot="plsmodel", sim="ordinary", stype="i", index=1)
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-family", family=binomial), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), typeboot="plsmodel", sim="ordinary", stype="i", index=1)

# Bastien CSDA 2005 Bootstrap

aze_compl.tilt.boot2 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic", family=NULL), statistic=coefs.plsRnp, R=c(499, 100, 100), alpha=c(0.025, 0.975), typeboot="fmodel_np", sim="ordinary", stype="i", index=1)
aze_compl.tilt.boot2 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic"), statistic=coefs.plsRnp, R=c(499, 100, 100), alpha=c(0.025, 0.975), typeboot="fmodel_np", sim="ordinary", stype="i", index=1)
aze_compl.tilt.boot2 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-family", family=binomial), statistic=coefs.plsRnp, R=c(499, 100, 100), alpha=c(0.025, 0.975), typeboot="fmodel_np", sim="ordinary", stype="i", index=1)



