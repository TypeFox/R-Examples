PseudoR2<-function(glmModel){
#if (length(class(object)) > 1)
#if (class(glmModel)!="glm" | class(glmModel)!="lm"){
#stop("Object not of class 'glm'")
#}
#else {
#	glmModel<-glmModel
#	}

logLikN<-glmModel$null/-2  ##log likelihood, null model
logLikF<-glmModel$dev/-2  ##log likelihood, full model
G2<-glmModel$null - glmModel$deviance
n <- length(glmModel$y)	
ystar <- predict(glmModel, type="response") 
class1<-ifelse(ystar >.5,1,0) 
classtab<-table(class1, glmModel$y, dnn=c("Predicted", "Actual")) ; maxOut<-max(margin.table(classtab, 2))
p<-glmModel$rank
penaltyN<- 2*(p)*(p+1)  ; penaltyD<- n-p-1  ; penalty<-penaltyN/penaltyD
MZystar <- predict(glmModel); sse <- sum((MZystar - mean(MZystar))^2) ; s2 <- switch(glmModel$family$link, "probit" = 1, "logit" = pi^2/3, NA)  #Needed for MZ R2
Enum<-sum((glmModel$y - ystar)^2); Edenom<-sum((glmModel$y - mean(glmModel$y))^2) #Needed for Effron R2

#R2s
r2McF<-1-logLikF/logLikN  #Mcfadden's R2
r2McFA<-1-(logLikF - p-1 )/logLikN #Mcfadden's Adj R2
r2CS<-1-exp(-G2/n) #ML Cox/Snell R2
r2N<-(1 - exp((glmModel$dev - glmModel$null)/n))/(1 - exp(-glmModel$null/n))# Nagelkerke/Cragg-Uhler R2
r2MZ<-sse / (n * s2 + sse)  #McKelvey and Zavoina pseudo R^2, using either the logit or probit link
r2E<-1-(Enum/Edenom) #Effron R2
r2C<-(classtab[1] + classtab[4])/n##Count R2 (proportion correctly classified)
r2CA<-(classtab[1] + classtab[4] - maxOut)/(n - maxOut) ##Adjusted Count R2 (proportion correctly classified)
aic<-2*(p)+glmModel$dev # AIC
Caic<-aic + penalty # AIC with a correction for finite sample size; useful with small sample sizes or a lot of predictors
 
results<-c(McFadden=r2McF, Adj.McFadden=r2McFA, Cox.Snell=r2CS, Nagelkerke=r2N, McKelvey.Zavoina=r2MZ, Effron=r2E, Count=r2C, Adj.Count=r2CA, AIC=aic, Corrected.AIC=Caic)
return(results)

}
