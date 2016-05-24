### R code from vignette source 'modTempEff.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("modTempEff")
options(prompt = "R> ", continue = "+   ")


###################################################
### code chunk number 2: loadlib
###################################################
library("modTempEff")


###################################################
### code chunk number 3: data
###################################################
data("dataDeathTemp")
head(dataDeathTemp)


###################################################
### code chunk number 4: modTempEff.Rnw:323-326 (eval = FALSE)
###################################################
## layout(matrix(c(1,1,2), ncol = 3))
## with(dataDeathTemp, plot(dec1, xlab="day", ylab="no. of deaths"))
## with(dataDeathTemp, plot(mtemp, dec1, xlab="temperature", ylab="no. of deaths"))


###################################################
### code chunk number 5: descrPlots
###################################################
layout(matrix(c(1,1,2), ncol = 3))
with(dataDeathTemp, plot(dec1, xlab="day", ylab="no. of deaths"))
with(dataDeathTemp, plot(mtemp, dec1, xlab="temperature", ylab="no. of deaths"))


###################################################
### code chunk number 6: modTempEff.Rnw:352-354
###################################################
o <- tempeff(dec1 ~ day + factor(dweek) + factor(year) + factor(month) + 
csdl(mtemp, psi = 20, L = c(60, 60)), data = dataDeathTemp, fcontrol = fit.control(display = TRUE))


###################################################
### code chunk number 7: modTempEff.Rnw:362-363
###################################################
o.noRidge <- update(o, .~. - day - factor(year) - factor(month) + seas(day, 30), fcontrol = fit.control(display = FALSE))


###################################################
### code chunk number 8: modTempEff.Rnw:377-380
###################################################
o

o.noRidge


###################################################
### code chunk number 9: modTempEff.Rnw:400-401
###################################################
o.Ridge.l <- update(o.noRidge, ridge = list(cold = "l", heat = "l"))


###################################################
### code chunk number 10: modTempEff.Rnw:404-405
###################################################
formula(o.Ridge.l)


###################################################
### code chunk number 11: modTempEff.Rnw:420-422
###################################################
o.Ridge.l4 <- update(o.noRidge, ridge = list(cold = "l^4", heat = "l^4"))
o.Ridge.l4


###################################################
### code chunk number 12: modTempEff.Rnw:428-429
###################################################
anova(o.noRidge, o.Ridge.l, o.Ridge.l4, test = "Cp")


###################################################
### code chunk number 13: modTempEff.Rnw:441-442
###################################################
summary(o.Ridge.l4)


###################################################
### code chunk number 14: modTempEff.Rnw:467-468
###################################################
coef(o.Ridge.l4,L=7)


###################################################
### code chunk number 15: DLcurve
###################################################
par(mfcol = c(2, 3))
plot(o.noRidge, new = FALSE)
plot(o.Ridge.l, new = FALSE)
plot(o.Ridge.l4, new = FALSE)


###################################################
### code chunk number 16: modTempEff.Rnw:504-508 (eval = FALSE)
###################################################
## par(mfcol = c(2, 3))
## plot(o.noRidge, new = FALSE)
## plot(o.Ridge.l, new = FALSE)
## plot(o.Ridge.l4, new = FALSE)


###################################################
### code chunk number 17: modTempEff.Rnw:523-525
###################################################
o2<-tempeff(dec1 ~ day + factor(dweek) + factor(year) + factor(month) + csdl(mtemp, psi = c(10, 20), L = c(60, 60)), data = dataDeathTemp,
     fcontrol = fit.control(display = TRUE))


###################################################
### code chunk number 18: modTempEff.Rnw:536-537
###################################################
summary(o2)


###################################################
### code chunk number 19: modTempEff.Rnw:545-548
###################################################
o2.Ridge.l4 <- update(o.Ridge.l4, psi = c(10, 20), fcontrol = fit.control(it.max = 30))

o2.Ridge.l4


###################################################
### code chunk number 20: modTempEff.Rnw:556-557
###################################################
anova(o.Ridge.l4, o2.Ridge.l4, test="BIC")


