### R code from vignette source 'Ch9.rnw'

###################################################
### code chunk number 1: setup
###################################################
source("GenericSettings.R")


###################################################
### code chunk number 2: Ch9.rnw:24-26
###################################################
data(HandHeight, package="MindOnStats")
str(HandHeight)


###################################################
### code chunk number 3: 3HandsVHeights
###################################################
plot(HandSpan~Height, data=HandHeight, ylab="Hand span")


###################################################
### code chunk number 4: Ch9.rnw:39-40
###################################################
summary(Hand.lm <- lm(HandSpan~Height, data=HandHeight))


###################################################
### code chunk number 5: 4HandResidsQQNorm
###################################################
Residuals = residuals(Hand.lm)
qqnorm(Residuals)
qqline(Residuals)


###################################################
### code chunk number 6: 4HandResidsVFits
###################################################
Fits=fitted(Hand.lm)
plot(Residuals~Fits)


###################################################
### code chunk number 7: XHandHeightFittedLinePlot
###################################################
plot(HandSpan~Height, data=HandHeight)
abline(Hand.lm)


###################################################
### code chunk number 8: Ch9.rnw:71-74
###################################################
data(Cereals, package="MindOnStats")
str(Cereals)
summary(Energy.lm <- lm(Energy~Dfibre, data=Cereals))


###################################################
### code chunk number 9: 6EnergyResidsQQNorm
###################################################
qqnorm(residuals(Energy.lm))
qqline(residuals(Energy.lm))


###################################################
### code chunk number 10: 6EnergyResidsVFits
###################################################
plot(residuals(Energy.lm)~fitted(Energy.lm))


###################################################
### code chunk number 11: Ch9.rnw:100-102
###################################################
data(Textbooks, package="MindOnStats")
summary(BookWeight.lm <- lm(Weight~Thickness, data=Textbooks))


###################################################
### code chunk number 12: 13BookWeightResidsQQNorm
###################################################
Residuals = residuals(BookWeight.lm)
qqnorm(Residuals)
qqline(Residuals)


###################################################
### code chunk number 13: 13BookWeightResidsVFits
###################################################
Fits = fitted(BookWeight.lm)
plot(Residuals~Fits)
abline(h=0)


###################################################
### code chunk number 14: Ch9.rnw:127-128
###################################################
sqrtWeight = sqrt(Textbooks$Weight)


###################################################
### code chunk number 15: Ch9.rnw:136-137
###################################################
summary(Energy.lm2 <- lm(Energy~Protein+Carbo+Dfibre, data=Cereals))


###################################################
### code chunk number 16: 15Energy4ResidGraphs
###################################################
Residuals=residuals(Energy.lm2)
Fits=fitted(Energy.lm2)
par(mfrow=c(2,2))
plot(Residuals~Cereals$Protein, xlab="Protein")
abline(h=0)
plot(Residuals~Cereals$Carbo, xlab="Carbohydrates")
abline(h=0)
plot(Residuals~Cereals$Dfibre, xlab="Dietary fibre")
abline(h=0)
plot(Residuals~Fits, xlab="Fitted values")
abline(h=0)


###################################################
### code chunk number 17: Ch9.rnw:159-161
###################################################
data(TimePerception, package="MindOnStats")
summary(TenSec.lm <- lm(TenSec~FiveSec*Gender, data=TimePerception))


###################################################
### code chunk number 18: XFitsVX
###################################################
plot(fitted(TenSec.lm)~TimePerception$FiveSec, xlab="Five seconds", main = "Predicted values for Ten Seconds vs actual data for Five Seconds", ylab="Fitted values")


###################################################
### code chunk number 19: Ch9.rnw:175-176
###################################################
summary( Price.lm1 <- lm(Price~Thickness+Weight+Year+Coverstyle+Colour+CD, data=Textbooks))


###################################################
### code chunk number 20: 16TextbooksPriceResidsWeight
###################################################
Residuals=residuals(Price.lm1)
attach(Textbooks)
plot(Residuals~Weight)
abline(h=0)


###################################################
### code chunk number 21: 16TextbooksPriceResidsColour
###################################################
plot(Residuals~Colour)


###################################################
### code chunk number 22: 16TextbooksPriceResidsYear
###################################################
plot(Residuals~Year)
abline(h=0)
detach(Textbooks)


###################################################
### code chunk number 23: Ch9.rnw:204-205
###################################################
summary(Price.lm2 <- lm(Price~poly(Thickness,2,raw=TRUE)+poly(Weight,2,raw=TRUE)+Year+Coverstyle+Colour+CD, data=Textbooks))


###################################################
### code chunk number 24: Ch9.rnw:214-215
###################################################
anova(Energy.lm)


###################################################
### code chunk number 25: Ch9.rnw:219-220
###################################################
confint(Energy.lm)


###################################################
### code chunk number 26: Ch9.rnw:227-229
###################################################
predict(Energy.lm, Cereals[1:5,], interval="prediction")
predict(Energy.lm, Cereals[1:5,], interval="confidence")


