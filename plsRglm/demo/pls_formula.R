set.seed(12345)

data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
PLS_glm_formula(log(ypine)~as.matrix(Xpine),3)
PLS_glm_formula(log(x11)~.,data=pine,3)
PLS_glm(dataY=log(ypine),dataX=Xpine,nt=3)

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=10,NK=1,modele="pls")
bbbbis <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,NK=1,modele="pls")
kfolds2CVinfos_glm(bbb)
kfolds2CVinfos_glm(bbbbis)

PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-gaussian",K=12)
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=TRUE,keepfolds=TRUE)$results_kfolds

#Different ways of model specifications
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian,K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds


bbb2 <- PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=10,modele="pls-glm-gaussian",keepcoeffs=TRUE)
bbb2 <- PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),K=6,keepcoeffs=TRUE)

#For Jackknife computations
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,1])

kfolds2Chisqind(bbb2)
kfolds2Chisq(bbb2)
kfolds2CVinfos_glm(bbb2)
PLS_lm(log(yCornell),XCornell,10,typeVC="standard")$InfCrit
rm(list=c("XCornell","yCornell","bbb","bbbbis","bbb2"))

data(Cornell)
PLS_lm_formula(Y~.,data=Cornell,10)$InfCrit
PLS_lm_formula(Y~.,data=Cornell,10,typeVC="standard")$InfCrit
PLS_lm_formula(Y~.,data=Cornell,6)$AIC 
PLS_lm_formula(Y~.,data=Cornell,6)$AIC.std    

data(Cornell)
PLS_lm_kfoldcv_formula(Y~.,data=Cornell,nt=3,K=12,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,data=Cornell,nt=3,K=12,keepfolds=FALSE)
PLS_lm_kfoldcv_formula(Y~.,data=Cornell,nt=3,K=6,NK=2,random=FALSE,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,data=Cornell,nt=3,K=6,NK=2,random=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,data=Cornell,nt=3,keepcoeffs=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,data=Cornell,nt=3,keepcoeffs=TRUE,keepfolds=FALSE)

bbb <- PLS_lm_kfoldcv_formula(Y~.,data=data.frame(scale(as.matrix(Cornell))[,]),nt=6,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv_formula(Y~.,data=data.frame(scale(as.matrix(Cornell))[,]),nt=6,K=6,NK=1)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb2)
PLS_lm_formula(Y~.,data=Cornell,6,typeVC="standard")$InfCrit
rm(list=c("bbb","bbb2"))

data(aze_compl)
# Lazraq-Cleroux PLS bootstrap Classic
aze_compl.boot <- bootplsglm(plsRglm(y~.,aze_compl,3,modele="pls-glm-logistic"), typeboot="plsmodel", sim="ordinary", stype="i", R=250)
boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))

# Bastien CSDA 2005 Bootstrap
aze_compl.boot2 <- bootplsglm(plsRglm(y~.,aze_compl,3,modele="pls-glm-logistic"), typeboot="fmodel_np", sim="ordinary", stype="i", R=250)
boxplots.bootpls(aze_compl.boot2)
confints.bootpls(aze_compl.boot2)
plots.confints.bootpls(confints.bootpls(aze_compl.boot2))
rm(list=c("aze_compl.boot","aze_compl.boot2"))

if("plsRcox" %in% installed.packages()){
# Only to provide an example of use of a factor - should be modeled with survival models
library(plsRcox)
data(Xmicro.censure_compl_imp)
data(micro.censure)
micro.censure.factor <- cbind(as.data.frame(Xmicro.censure_compl_imp[,-40]),STADE=factor(Xmicro.censure_compl_imp[,40]),survyear=micro.censure$survyear)
str(micro.censure.factor)
PLS_lm_formula(survyear~.,micro.censure.factor)$Coeffs

# Lazraq-Cleroux PLS bootstrap Classic
micro.boot <- bootpls(plsR(survyear~STADE,micro.censure.factor,nt=3,modele="pls"), typeboot="plsmodel", sim="ordinary", stype="i", R=250)
boxplots.bootpls(micro.boot,indices=2:6)
# Confidence intervals plotting
confints.bootpls(micro.boot,indices=2:6)
plots.confints.bootpls(confints.bootpls(micro.boot,indices=2:6))

# Bastien CSDA 2005 Bootstrap
micro.boot2 <- bootpls(plsR(survyear~STADE,micro.censure.factor,nt=3,modele="pls"), typeboot="fmodel_np", sim="ordinary", stype="i", R=250)
boxplots.bootpls(micro.boot2,indices=2:6)
# Confidence intervals plotting
confints.bootpls(micro.boot2,indices=2:6)
plots.confints.bootpls(confints.bootpls(micro.boot2,indices=2:6))
}






