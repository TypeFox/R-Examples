set.seed(12345)

data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
modpls <- PLS_glm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit
modpls <- plsRglm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit


XbordeauxNA<-Xbordeaux
XbordeauxNA[1,1] <- NA
modplsNA <- PLS_glm(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit
modplsNA <- plsRglm(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
dataset <- cbind(y=ybordeaux,Xbordeaux)
options(contrasts = c("contr.treatment", "contr.poly"))
bordeaux.boot3<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="permutation", stype="i", R=250)
boxplots.bootpls(bordeaux.boot3)
boxplots.bootpls(bordeaux.boot3,ranget0=TRUE)

bordeaux.boot2<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="permutation", stype="i", R=250, strata=unclass(ybordeaux))
boxplots.bootpls(bordeaux.boot2)


bordeaux.boot4<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="balanced", stype="i", R=250)
boxplots.bootpls(bordeaux.boot4)
boot.ci(bordeaux.boot4)
confints.bootpls(bordeaux.boot4)
plots.confints.bootpls(confints.bootpls(bordeaux.boot4))

bordeaux.boot5<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="balanced", stype="i", R=250, strata=unclass(ybordeaux))
boxplots.bootpls(bordeaux.boot5)
#boot.ci(bordeaux.boot5)
confints.bootpls(bordeaux.boot5)
plots.confints.bootpls(confints.bootpls(bordeaux.boot5))


if(require(chemometrics)){
data(hyptis)
hyptis
yhyptis <- factor(hyptis$Group,ordered=TRUE)
Xhyptis <- as.data.frame(hyptis[,c(1:6)])
dataset <- cbind(y=yhyptis,Xhyptis)
options(contrasts = c("contr.treatment", "contr.poly"))
modpls2 <- plsRglm(yhyptis,Xhyptis,6,modele="pls-glm-polr")
modpls2$Coeffsmodel_vals
modpls2$InfCrit
modpls2$Coeffs
modpls2$std.coeffs

table(yhyptis,predict(modpls2$FinalModel,type="class"))


modpls3 <- PLS_glm(yhyptis[-c(1,2,3)],Xhyptis[-c(1,2,3),],3,modele="pls-glm-polr",dataPredictY=Xhyptis[c(1,2,3),])
bbb <- PLS_glm_kfoldcv(yhyptis,Xhyptis,nt=4,K=10,random=TRUE,modele="pls-glm-polr",keepcoeffs=T)
kfolds2CVinfos_glm(bbb,MClassed=TRUE)


# Lazraq-Cleroux PLS bootstrap Classic
hyptis.boot0<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), typeboot="plsmodel", sim="ordinary", stype="i", R=250)
rownames(hyptis.boot0$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot0)
boxplots.bootpls(hyptis.boot0,xaxisticks=FALSE)
boxplots.bootpls(hyptis.boot0,ranget0=TRUE)
boxplots.bootpls(hyptis.boot0,ranget0=TRUE,xaxisticks=FALSE)


hyptis.boot1<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), typeboot="plsmodel", sim="permutation", stype="i", R=250)
rownames(hyptis.boot1$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot1)
boxplots.bootpls(hyptis.boot1,xaxisticks=FALSE)
boxplots.bootpls(hyptis.boot1,ranget0=TRUE)
boxplots.bootpls(hyptis.boot1,ranget0=TRUE,xaxisticks=FALSE)


hyptis.boot2<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), typeboot="plsmodel", sim="balanced", stype="i", R=250)
rownames(hyptis.boot2$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot2)
boxplots.bootpls(hyptis.boot2,xaxisticks=FALSE)
confints.bootpls(hyptis.boot2)
plots.confints.bootpls(confints.bootpls(hyptis.boot2),legendpos = "bottomleft",xaxisticks=TRUE)
plots.confints.bootpls(confints.bootpls(hyptis.boot2),legendpos = "bottomleft",xaxisticks=FALSE)


# Bastien CSDA 2005 Bootstrap
hyptis.boot3<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), typeboot="fmodel_np", sim="ordinary", stype="i", R=250)
rownames(hyptis.boot3$t0)<-c("Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot3)
boxplots.bootpls(hyptis.boot3,xaxisticks=FALSE)
boxplots.bootpls(hyptis.boot3,ranget0=TRUE)
boxplots.bootpls(hyptis.boot3,ranget0=TRUE,xaxisticks=FALSE)


hyptis.boot4<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), typeboot="fmodel_np", sim="permutation", stype="i", R=250)
rownames(hyptis.boot4$t0)<-c("Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot4)
boxplots.bootpls(hyptis.boot4,xaxisticks=FALSE)
boxplots.bootpls(hyptis.boot4,ranget0=TRUE)
boxplots.bootpls(hyptis.boot4,ranget0=TRUE,xaxisticks=FALSE)


hyptis.boot5<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), typeboot="fmodel_np", sim="balanced", stype="i", R=250)
rownames(hyptis.boot5$t0)<-c("Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot5)
boxplots.bootpls(hyptis.boot5,xaxisticks=FALSE)
confints.bootpls(hyptis.boot5)
plots.confints.bootpls(confints.bootpls(hyptis.boot5),legendpos = "bottomleft",xaxisticks=TRUE)
plots.confints.bootpls(confints.bootpls(hyptis.boot5),legendpos = "bottomleft",xaxisticks=FALSE)
}
