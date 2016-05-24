################################################################################
##                                                                            ##
##                               mvMORPH: test.LRT                            ##
##                                                                            ##
##  Created by Julien Clavel - 02-02-2015                                     ##
##  (julien.clavel@hotmail.fr/ clavel@biologie.ens.fr)                        ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################

LRT<-function(model1,model2,echo=TRUE){
##-------------------LRT comparison of the models-----------------------------##

if(any(class(model1)[2]!=class(model2)[2])){warning("You are using an LRT test on non-nested models!! Use AIC comparison instead")}

if(model1$param$nparam<=model2$param$nparam){
    mod1<-model1;mod2<-model2
    L2<-model1$LogLik
    L1<-model2$LogLik
    model1<-mod2;model2<-mod1
    warning("The order of the models has been changed for the comparison")
}else{
L1<-model1$LogLik
L2<-model2$LogLik
}
# Set names
if(class(model1)[2]=="bm" | class(model2)[2]=="bm"){
    model1$param$model<-if(model1$param$constraint==TRUE){paste(model1$param$model," constrained")}else if(model1$param$constraint=="shared"){paste(model1$param$model," shared variances")}else if(model1$param$constraint=="proportional"){paste(model1$param$model," proportional")}else if(model1$param$constraint=="correlation"){paste(model1$param$model," equal correlation")}else if(model1$param$constraint=="diagonal"){paste(model1$param$model," diagonal")}else{model1$param$model}

    model2$param$model<-if(model2$param$constraint==TRUE){paste(model2$param$model," constrained")}else if(model2$param$constraint=="shared"){paste(model2$param$model," shared variances")}else if(model2$param$constraint=="proportional"){paste(model2$param$model," proportional")}else if(model2$param$constraint=="correlation"){paste(model2$param$model," equal correlation")}else if(model1$param$constraint=="diagonal"){paste(model1$param$model," diagonal")}else{model2$param$model}
}

if(class(model1)[2]=="ou" | class(model2)[2]=="ou"){
    model1$param$model<-if(model1$param$decomp=="diagonal"){paste(model1$param$model," diagonal")}else if(model1$param$decomp=="symmetricPositive"){paste(model1$param$model," symmetric positive")}else if(model1$param$decomp=="symmetric"){paste(model1$param$model," symmetric")}else if(model1$param$decomp=="nsymmetric"){paste(model1$param$model," non-symmetric")}else if(model1$param$decomp=="nsymPositive"){paste(model1$param$model," non-symmetric positive")}else{model1$param$model}
        
    model2$param$model<-if(model2$param$decomp=="diagonal"){paste(model2$param$model," diagonal")}else if(model2$param$decomp=="symmetricPositive"){paste(model2$param$model," symmetric positive")}else if(model2$param$decomp=="symmetric"){paste(model2$param$model," symmetric")}else if(model2$param$decomp=="nsymmetric"){paste(model2$param$model," non-symmetric")}else if(model2$param$decomp=="nsymPositive"){paste(model2$param$model," non-symmetric positive")}else{model2$param$model}
}

#
LRT<-(2*((L1-L2)))
#difference in degrees of freedom
ddf<-model1$param$nparam-model2$param$nparam
LRT.prob<-pchisq(LRT,ddf,lower.tail=FALSE)
if(echo==TRUE){
    if(LRT.prob<0.001){signif<-c("***")}else if(LRT.prob<0.01){
        signif<-c("**") }else if(LRT.prob<0.05){signif<-c("*")}else if(LRT.prob<0.1){signif<-c(".")}else{signif<-""}
    cat("-- Log-likelihood Ratio Test --","\n")
    cat("Model",model1$param$model," versus ",model2$param$model,"\n")
    cat("Number of degrees of freedom :",ddf,"\n")
    cat("LRT statistic:",LRT," p-value:",LRT.prob,signif,"\n")
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}

results<-list(pval=LRT.prob, ratio=LRT, ddf=ddf, model1=model1$param$model, model2=model2$param$model)
class(results)<-c("mvmorph.lrt")
invisible(results)
#End
}
