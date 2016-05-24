### R code from vignette source 'kf.Rnw'

###################################################
### code chunk number 1: kf.Rnw:258-259
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: kf.Rnw:261-262 (eval = FALSE)
###################################################
## install.packages("cts")


###################################################
### code chunk number 3: kf.Rnw:265-267 (eval = FALSE)
###################################################
## library("cts")
## vignette("kf",package = "cts")


###################################################
### code chunk number 4: kf.Rnw:270-271 (eval = FALSE)
###################################################
## edit(vignette("kf",package = "cts"))


###################################################
### code chunk number 5: kf.Rnw:279-281
###################################################
library("cts")
data("V22174")


###################################################
### code chunk number 6: kf.Rnw:287-290
###################################################
### plot the oxygen isotope time series
plot(V22174,type="l",xlab="Time in kiloyears", ylab="")
rug(V22174[,1], col="red")


###################################################
### code chunk number 7: kf.Rnw:295-296
###################################################
time <- system.time(V22174.car14 <- car(V22174,scale=0.2,order=14))[1]


###################################################
### code chunk number 8: kf.Rnw:300-304
###################################################
### fit the modified form of a continuous AR(14) model
V22174.car14 <- car(V22174,scale=0.2,order=14)
tab1 <- cbind(V22174.car14$tnit, V22174.car14$ss, V22174.car14$bit[,14])
colnames(tab1) <- c("Iteration","Sum of Squares","phi_14")


###################################################
### code chunk number 9: kf.Rnw:307-308
###################################################
print(as.data.frame(round(tab1,5)),row.names=FALSE, print.gap=8)


###################################################
### code chunk number 10: kf.Rnw:312-314
###################################################
### AIC output based on Belcher et. al (1994)
AIC(V22174.car14)


###################################################
### code chunk number 11: kf.Rnw:316-319
###################################################
### fit the modified form of a continuous AR(7) model
V22174.car7 <- car(V22174,scale=0.2,order=7)
summary(V22174.car7)


###################################################
### code chunk number 12: kf.Rnw:323-333
###################################################
### classical AIC and BIC results for model selection
norder <- 14
V22174.aic <- V22174.bic <- rep(NA, norder)
for (i in 1:norder){
fit <- car(V22174,scale=0.2,order=i)
V22174.aic[i] <- fit$aic
V22174.bic[i] <- fit$bic
}
res <- data.frame(order=1:norder, AIC=V22174.aic, BIC=V22174.bic)
print(res, row.names=FALSE, print.gap=8)


###################################################
### code chunk number 13: kf.Rnw:339-343
###################################################
### plot the spectrum for the modified contiuous AR(14) and AR(7) models
par(mfrow=c(2,1))
spectrum(V22174.car14)
spectrum(V22174.car7)


###################################################
### code chunk number 14: kf.Rnw:361-363
###################################################
### model diagnostics check
tsdiag(V22174.car7)


###################################################
### code chunk number 15: kf.Rnw:372-373
###################################################
data("asth")


###################################################
### code chunk number 16: kf.Rnw:375-377 (eval = FALSE)
###################################################
## plot(asth,type="l",xlab="Time in hours", ylab="")
## rug(asth[,1], col="red")


###################################################
### code chunk number 17: kf.Rnw:382-385
###################################################
### plot of the lung function time series
plot(asth,type="l",xlab="Time in hours", ylab="")
rug(asth[,1], col="red")


###################################################
### code chunk number 18: kf.Rnw:391-394
###################################################
### fit the modified form of a continuous AR(4) model
asth.car4 <- car(asth,scale=0.25,order=4, ctrl=car_control(n.ahead=10))
summary(asth.car4)


###################################################
### code chunk number 19: kf.Rnw:401-404
###################################################
### fit the modified form of a continuous AR(7) model with measurement error
asth.vri <- car(asth,scale=0.25,order=4, ctrl=car_control(vri=TRUE))
summary(asth.vri)


###################################################
### code chunk number 20: kf.Rnw:409-413
###################################################
### plot of spectrum for the continuous AR(4) model without/with measurement error
par(mfrow=c(2,1))
spectrum(asth.car4)
spectrum(asth.vri)


###################################################
### code chunk number 21: kf.Rnw:420-422
###################################################
### determine the zeros of equation (3)
factab(asth.car4)


###################################################
### code chunk number 22: kf.Rnw:428-435
###################################################
### decompose the original time series into three corresponding components 
### via the Kalman smoother 
asth.kalsmo <- kalsmo(asth.car4)
par(mfrow=c(3,1))
kalsmoComp(asth.kalsmo,comp=1, xlab="Time in hours")
kalsmoComp(asth.kalsmo,comp=c(2,3), xlab="Time in hours")
kalsmoComp(asth.kalsmo,comp=4,xlab="Time in hours")


###################################################
### code chunk number 23: kf.Rnw:445-447
###################################################
### predict the last 10 steps past the end of time series
predict(asth.car4, xlab="Time in hours")


