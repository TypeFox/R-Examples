### R code from vignette source 'penalized.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(continue = "  ")


###################################################
### code chunk number 2: load
###################################################
library(penalized)
library(survival)
data(nki70)


###################################################
### code chunk number 3: setseed
###################################################
set.seed(1)


###################################################
### code chunk number 4: first
###################################################
fit <- penalized(ER, ~DIAPH3+NUSAP1, data=nki70, lambda2=1)
fit <- penalized(ER, nki70[,10:11], data=nki70, lambda2=1)
fit <- penalized(ER~DIAPH3+NUSAP1, data=nki70, lambda2=1)


###################################################
### code chunk number 5: survival
###################################################
fit <- penalized(Surv(time,event)~DIAPH3+NUSAP1, data=nki70, lambda2=1)


###################################################
### code chunk number 6: attach
###################################################
attach(nki70)


###################################################
### code chunk number 7: elastic net
###################################################
fit <- penalized(Surv(time,event)~DIAPH3+NUSAP1, data=nki70, lambda1=1, lambda2=1)


###################################################
### code chunk number 8: vector lambda
###################################################
fit <- penalized(Surv(time,event)~DIAPH3+NUSAP1, data=nki70, lambda2=c(1,2))


###################################################
### code chunk number 9: extract
###################################################
residuals(fit)[1:10]
fitted(fit)[1:10]
basesurv(fit)


###################################################
### code chunk number 10: coefficients
###################################################
coefficients(fit, "all")


###################################################
### code chunk number 11: loglik_penalty
###################################################
loglik(fit)
penalty(fit)


###################################################
### code chunk number 12: predict
###################################################
predict(fit, ~DIAPH3+NUSAP1, data=nki70[1:3,])
predict(fit, nki70[1:3,c("DIAPH3","NUSAP1")])


###################################################
### code chunk number 13: predict_survival
###################################################
pred <- predict(fit, nki70[1:3,c("DIAPH3","NUSAP1")])
survival(pred, time=5)


###################################################
### code chunk number 14: weights
###################################################
coefficients(fit)
coefficients(fit, standardize = TRUE)
weights(fit)


###################################################
### code chunk number 15: unpenalized
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], ~ER, lambda2=1)
fit <- penalized(Surv(time,event)~ER, nki70[,8:77], lambda2=1)


###################################################
### code chunk number 16: strata
###################################################
fit <- penalized(Surv(time,event)~strata(ER), nki70[,8:77], lambda2=1)


###################################################
### code chunk number 17: steps
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], lambda1=1,
    steps=50, trace = FALSE)
plotpath(fit, log="x")


###################################################
### code chunk number 18: steps park
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], lambda1=1,
    steps="Park", trace = FALSE)


###################################################
### code chunk number 19: stepsplot
###################################################
plotpath(fit, log="x")


###################################################
### code chunk number 20: positive
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], positive=TRUE)


###################################################
### code chunk number 21: positivestepsplot
###################################################
fit0 <- penalized(Surv(time,event), nki70[,8:77], positive=TRUE,
    steps=50)
plotpath(fit0)


###################################################
### code chunk number 22: positivecoefficients
###################################################
coefficients(fit)


###################################################
### code chunk number 23: partpositive
###################################################
coef(penalized(Surv(time,event), nki70[,8:16], positive=c(F,rep(T,8))))


###################################################
### code chunk number 24: fused lasso
###################################################
 X <- matrix(0,70,100)
  for (i in 1:100){
     X[1:70,i] <- rnorm(70,mean=0,sd=0.5)
   }
 colnames(X) = as.character (1:ncol(X))
  rownames(X) = as.character (1:nrow(X))


###################################################
### code chunk number 25: fused lasso
###################################################
 a <- sample(1:ncol(X),50,prob=rep(0.5,length(1:ncol(X))))
 for (i in 1:50){
      X[30:40,a[i]]<-rnorm(length(30:40),mean = -0.7 ,sd=0.5)
    }


###################################################
### code chunk number 26: fused lasso
###################################################
    Xbeta <- rnorm(100, mean = 0, sd = 0.5)
    Xbeta[a] <- rnorm (length(a) , mean = -0.7 , sd = 0.5)

    beta <- numeric(100)
    beta [1:length(beta)] <- -1

   responsep <- numeric(100)

 for(i in 1:100){
     coeff <- -beta[i] * Xbeta[i]
     responsep[i] <- 1/(1+exp(coeff))
   }

  X <- t(X)
    response=responsep
      for(i in 1:100){
      response[i] <- sample(0:1 , size = 1  , prob = c((1-responsep[i]),responsep[i]))
  }


###################################################
### code chunk number 27: fused lasso
###################################################
fit <- penalized(response, X, lambda1 = 2, lambda2=2,fusedl=TRUE)


###################################################
### code chunk number 28: fusedlassocoeffplot
###################################################
fit <- penalized(response, X, lambda1 = 2, lambda2=3,fusedl=TRUE)
plot(coefficients(fit,"all")[-1],main = "fused lasso", col="red",xlab = "probes",ylab = "coefficients",type="l")


###################################################
### code chunk number 29: fused lasso
###################################################
chr = c(rep(1,30),rep(2,20),rep(3,10),rep(4,10))
fit <- penalized(response, X, lambda1 = 2, lambda2=2,fusedl=chr)


###################################################
### code chunk number 30: globaltest_install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("globaltest")


###################################################
### code chunk number 31: globaltest
###################################################
library(globaltest)
gt(Surv(time,event), nki70[,8:77])


###################################################
### code chunk number 32: cvl1
###################################################
fit <- cvl(Surv(time,event), nki70[,8:77], lambda1=1, fold=10)


###################################################
### code chunk number 33: cvl2
###################################################
fit$cvl
fit$fullfit


###################################################
### code chunk number 34: cvl3
###################################################
fit <- cvl(Surv(time,event), nki70[,8:77], lambda1=2, fold=fit$fold)


###################################################
### code chunk number 35: cvlappr1
###################################################
fit1 <- cvl(Surv(time,event), nki70[,8:77], lambda2=10)


###################################################
### code chunk number 36: cvlappr2
###################################################
fit1$cvl


###################################################
### code chunk number 37: cvlappr3
###################################################
fit2 <- cvl(Surv(time,event), nki70[,8:77], lambda2=10, approximate=TRUE)


###################################################
### code chunk number 38: cvlappr4
###################################################
fit2$cvl


###################################################
### code chunk number 39: breslow
###################################################
fit$predictions
time(fit$predictions)
as.data.frame(basesurv(fit$fullfit))[1:10,]
plot(fit$predictions)


###################################################
### code chunk number 40: breslowplot
###################################################
plot(fit$predictions)


###################################################
### code chunk number 41: cv-survival
###################################################
survival(fit$predictions, 5)[1:10]


###################################################
### code chunk number 42: prof
###################################################
fit1 <- profL1(Surv(time,event), nki70[,50:70],fold=10, plot=TRUE)
fit2 <- profL2(Surv(time,event), nki70[,50:70],fold=fit1$fold,
    minl = 0.01, maxl = 1000)
plot(fit2$lambda, fit2$cvl, type="l", log="x")


###################################################
### code chunk number 43: profplot1
###################################################
plot(fit1$lambda, fit1$cvl, type="l")


###################################################
### code chunk number 44: profplot2
###################################################
plot(fit2$lambda, fit2$cvl, type="l", log="x")


###################################################
### code chunk number 45: profpath
###################################################
plotpath(fit2$fullfit, log="x")


###################################################
### code chunk number 46: profpathplot
###################################################
plotpath(fit2$fullfit, log="x")


###################################################
### code chunk number 47: opt1
###################################################
opt1 <- optL1(Surv(time,event), nki70[,50:70], fold=fit1$fold)


###################################################
### code chunk number 48: optres
###################################################
opt1$lambda
opt1$cvl


###################################################
### code chunk number 49: opt2
###################################################
opt2 <- optL2(Surv(time,event), nki70[,50:70], fold=fit2$fold)


