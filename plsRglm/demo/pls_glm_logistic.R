set.seed(12345)

data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
modpls$valpvalstep
modpls$Coeffsmodel_vals

plot(plsRglm(yaze_compl,Xaze_compl,4,modele="pls-glm-logistic")$FinalModel)
plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$pvalstep
plot(plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$FinalModel)
rm(list=c("Xaze_compl","yaze_compl","modpls"))




dimX <- 6
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 2
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar2)[,1]
Xsimbin1 <- dicho(dataAstar2)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar2","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 3
dataAstar3 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar3)[,1]
Xsimbin1 <- dicho(dataAstar3)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar3","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 5
dataAstar5 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar5)[,1]
Xsimbin1 <- dicho(dataAstar5)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar5","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 6
dataAstar6 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar6)[,1]
Xsimbin1 <- dicho(dataAstar6)[,2:(dimX+1)]
modplsglm <- plsRglm(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit




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

plot(PLS_glm(yaze_compl,Xaze_compl,4,modele="pls-glm-logistic")$FinalModel)
PLS_glm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$pvalstep
plot(PLS_glm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$FinalModel)



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


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
PLS_glm(yaze_compl,Xaze_compl,10,modele="pls-glm-logistic",typeVC="none")$InfCrit
PLS_glm_wvc(yaze_compl,Xaze_compl,10,modele="pls-glm-logistic", keepcoeffs=TRUE)



data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
plsRglmmodel.default(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit

bbb <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls")
kfolds2CVinfos_glm(bbb,MClassed=TRUE)
bbba <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls-glm-family",family=gaussian())
kfolds2CVinfos_glm(bbba,MClassed=TRUE)
bbb2 <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls-glm-logistic")
kfolds2CVinfos_glm(bbb2,MClassed=TRUE)
bbb2a <- PLS_glm_kfoldcv(yaze_compl,Xaze_compl,nt=10,K=12,NK=1,keepfolds=FALSE,keepdataY=TRUE,modele="pls-glm-family",family=binomial())
kfolds2CVinfos_glm(bbb2a,MClassed=TRUE)

bbb.lo <- PLS_glm_kfoldcv(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-family",family="binomial")
bbb2.lo <- PLS_glm_kfoldcv(dataY=yaze_compl,dataX=Xaze_compl,nt=4,modele="pls-glm-logistic")
kfolds2Chisq(bbb.lo)
kfolds2Chisq(bbb2.lo)


aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), typeboot="fmodel_np", sim="ordinary", stype="i", R=1000)

(temp.ci <- confints.bootpls(aze_compl.boot3,1:33,typeBCa=FALSE))
(temp.ci <- confints.bootpls(aze_compl.boot3,c(2,4,6)))
(temp.ci <- confints.bootpls(aze_compl.boot3))

boxplots.bootpls(aze_compl.boot3)
boxplots.bootpls(aze_compl.boot3,prednames=FALSE)
boxplots.bootpls(aze_compl.boot3,prednames=FALSE,articlestyle=FALSE,main="Bootstrap distribution for the bj")
boxplots.bootpls(aze_compl.boot3,indices=1:33,prednames=FALSE)
boxplots.bootpls(aze_compl.boot3,indices=c(2,4,6),prednames=FALSE)

plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,prednames=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,main="Bootstrap confidence intervals for the bj")
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)
plots.confints.bootpls(temp.ci,c(2,4,6),"bottomleft")
plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,main="Bootstrap confidence intervals for some of the bj")
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)

temp.ci <- confints.bootpls(aze_compl.boot3,1:33,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,indices=1:33,prednames=FALSE)



data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
dataset <- cbind(y=yaze_compl,Xaze_compl)
library(boot)
# Lazraq-Cleroux PLS bootstrap Classic
#aze_compl.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, nt=3, modele="pls-glm-logistic")
# The same
#aze_compl.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, nt=3, modele="pls-glm-family",family=binomial)
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), typeboot="plsmodel", sim="ordinary", stype="i", R=250)
# The same
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-family",family=binomial), typeboot="plsmodel", sim="ordinary", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)

boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))

plot(aze_compl.boot,index=2)
jack.after.boot(aze_compl.boot, index=2, useJ=TRUE, nt=3)
plot(aze_compl.boot, index=2,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic"), typeboot="plsmodel", statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)





# Bastien CSDA 2005 bootstrap
#aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), typeboot="fmodel_np", sim="ordinary", stype="i", R=250)
# The same
aze_compl.boot3 <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-family",family=binomial), typeboot="fmodel_np", sim="ordinary", stype="i", R=250)
boot.ci(aze_compl.boot3, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze_compl.boot3, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze_compl.boot3, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze_compl.boot3, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze_compl.boot3, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze_compl.boot3, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)

boxplots.bootpls(aze_compl.boot3)
confints.bootpls(aze_compl.boot3)
plots.confints.bootpls(confints.bootpls(aze_compl.boot3))

plot(aze_compl.boot3,index=2)
jack.after.boot(aze_compl.boot3, index=2, useJ=TRUE, nt=3)
plot(aze_compl.boot3, index=2,jack=TRUE)
aze_compl.tilt.boot2 <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic"), typeboot="fmodel_np", statistic=coefs.plsRnp, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)




















# PLS bootstrap balanced

aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), typeboot="fmodel_np", sim="balanced", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)


boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))



plot(aze_compl.boot)
jack.after.boot(aze_compl.boot, index=1, useJ=TRUE, nt=3)
plot(aze_compl.boot,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)


# PLS permutation bootstrap

#aze_compl.boot2 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", stype="i", R=250, nt=3, modele="pls-glm-logistic")
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), sim="permutation", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=6)
boxplots.bootpls(aze_compl.boot)
plot(aze_compl.boot)


#With missing data
data(aze)
Xaze<-aze[,2:34]
yaze<-aze$y

dataset <- cbind(y=yaze,Xaze)

library(boot)
#aze.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, nt=3, modele="pls-glm-logistic")
aze.boot <- bootplsglm(plsRglm(yaze,Xaze,3,modele="pls-glm-logistic"), sim="ordinary", stype="i", R=250)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)


boxplots.bootpls(aze.boot)
confints.bootpls(aze.boot)
plots.confints.bootpls(confints.bootpls(aze.boot))

plot(aze.boot)
jack.after.boot(aze.boot, index=1, useJ=TRUE, nt=3)
plot(aze.boot,jack=TRUE)


