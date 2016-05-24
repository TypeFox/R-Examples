### R code from vignette source 'HDPenReg.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: signal
###################################################
library(HDPenReg)
simu=simul(50,10000,0.4,2,750,matrix(c(0.1,0.9,0.0001,0.0001),nrow=2),10)
case=which(simu$response==1)
plot(simu$data[case[1],],ylim=c(0,4),pch=".",xlab="Value")


###################################################
### code chunk number 2: lars (eval = FALSE)
###################################################
## # HDlars(X, y, maxSteps,intercept,eps)


###################################################
### code chunk number 3: simdata1
###################################################
data1=simul(50,10000,0.4,2,1000,matrix(c(0,1,0,0),nrow=2),5)


###################################################
### code chunk number 4: exlars
###################################################
reslars=HDlars(data1$data, data1$response)


###################################################
### code chunk number 5: plotlars (eval = FALSE)
###################################################
## plot(reslars)


###################################################
### code chunk number 6: plotlarsb
###################################################
plot(reslars)


###################################################
### code chunk number 7: snpcaus
###################################################
sort(data1$causalSNP)
sort(reslars@variable[[16]])


###################################################
### code chunk number 8: cvlarsdef (eval = FALSE)
###################################################
## # HDcvlars(X, y, nbFolds, index, mode, maxSteps, partition, intercept, eps)


###################################################
### code chunk number 9: cvlars
###################################################
rescvlars=HDcvlars(data1$data, data1$response,5)


###################################################
### code chunk number 10: plotcvlars (eval = FALSE)
###################################################
## plot(rescvlars)


###################################################
### code chunk number 11: plotcvlarsb
###################################################
plot(rescvlars)


###################################################
### code chunk number 12: coeff
###################################################
coefficients=computeCoefficients(reslars, rescvlars$minIndex, mode = "lambda")


###################################################
### code chunk number 13: pred (eval = FALSE)
###################################################
## # yPred=predict(reslars, rescvlars$fraction, mode = "fraction")


###################################################
### code chunk number 14: fusion (eval = FALSE)
###################################################
## # HDfusion(X, y, maxSteps,intercept,eps)


###################################################
### code chunk number 15: exfusion
###################################################
resfusion=HDfusion(data1$data, data1$response)


###################################################
### code chunk number 16: plotfusion (eval = FALSE)
###################################################
## plot(resfusion)


###################################################
### code chunk number 17: plotfusion2
###################################################
plot(resfusion)


###################################################
### code chunk number 18: plotcoeff
###################################################
plotCoefficient(resfusion,20)


