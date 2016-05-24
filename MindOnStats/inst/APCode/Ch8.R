### R code from vignette source 'Ch8.rnw'

###################################################
### code chunk number 1: setup
###################################################
source("GenericSettings.R")


###################################################
### code chunk number 2: Ch8.rnw:21-23
###################################################
data(GoGoGo, package="MindOnStats")
str(GoGoGo)


###################################################
### code chunk number 3: Ch8.rnw:29-31
###################################################
GreenLights.aov = aov(Time~AgeGroup, data=GoGoGo, subset=Lights=="green")
summary(GreenLights.aov)


###################################################
### code chunk number 4: Ch8.rnw:36-37
###################################################
data(Textbooks, package="MindOnStats")


###################################################
### code chunk number 5: Ch8.rnw:41-42
###################################################
summary(Price.aov1 <- aov(Price~Discipline, data=Textbooks))


###################################################
### code chunk number 6: Ch8.rnw:51-54
###################################################
data(PaperPlanes, package="MindOnStats")
dim(PaperPlanes)
head(PaperPlanes)


###################################################
### code chunk number 7: Ch8.rnw:58-61
###################################################
StingRay = PaperPlanes[PaperPlanes$Design=="stingray glider",]
dim(StingRay)
head(StingRay)


###################################################
### code chunk number 8: 4StingRayBoxPlots
###################################################
boxplot(FlightTime~Paper, data=StingRay)


###################################################
### code chunk number 9: Ch8.rnw:73-74
###################################################
summary(StingRay.aov1 <- aov(FlightTime~Paper, data=StingRay))


###################################################
### code chunk number 10: XGreenLightsResidPlots
###################################################
par(mfrow=c(2,2))
plot(GreenLights.aov)


###################################################
### code chunk number 11: 6GreenLightsResidFits
###################################################
Residuals = residuals(GreenLights.aov)
Fits = fitted(GreenLights.aov)
plot(Residuals~Fits)


###################################################
### code chunk number 12: 6GreenLightsQQNorm
###################################################
qqnorm(Residuals)
qqline(Residuals)


###################################################
### code chunk number 13: Ch8.rnw:112-113
###################################################
summary(StingRay.aov2 <- aov(log(FlightTime)~Paper, data=StingRay))


###################################################
### code chunk number 14: 9StingRayResidQQNorm
###################################################
Residuals = residuals(StingRay.aov2)
qqnorm(Residuals)
qqline(Residuals)


###################################################
### code chunk number 15: 9StingRayResidVFits
###################################################
Fits = fitted(StingRay.aov2)
plot(Residuals~Fits)


###################################################
### code chunk number 16: 10StingRayLogResidQQNorm
###################################################
Residuals = residuals(StingRay.aov2)
qqnorm(Residuals)
qqline(Residuals)


###################################################
### code chunk number 17: 10StingRayLogResidVFits
###################################################
Fits = fitted(StingRay.aov2)
plot(Residuals~Fits)


###################################################
### code chunk number 18: Ch8.rnw:157-158
###################################################
kruskal.test(FlightTime~Paper, data=StingRay)


###################################################
### code chunk number 19: Ch8.rnw:170-172
###################################################
GreenLights.HSD = TukeyHSD(GreenLights.aov)
GreenLights.HSD


###################################################
### code chunk number 20: 13Tukey
###################################################
plot(GreenLights.HSD)


###################################################
### code chunk number 21: Ch8.rnw:188-189
###################################################
data(TimePerception, package="MindOnStats")


###################################################
### code chunk number 22: Ch8.rnw:193-194
###################################################
summary(aov(TenSec~Gender*AgeGroup, data=TimePerception))


###################################################
### code chunk number 23: 15TimeIntPlot
###################################################
attach(TimePerception)
interaction.plot(response=TenSec, x.factor=AgeGroup, trace.factor=Gender)
detach(TimePerception)


###################################################
### code chunk number 24: Ch8.rnw:210-213
###################################################
summary(Planes.aov2 <- aov(FlightTime~Paper*Design, data=PaperPlanes))
lnFlightTime = log(PaperPlanes$FlightTime)
summary(Planes.aov3 <- aov(lnFlightTime~Paper*Design, data=PaperPlanes))


###################################################
### code chunk number 25: 18PlanesResidQQNorm
###################################################
Residuals = residuals(Planes.aov2)
qqnorm(Residuals)
qqline(Residuals)


###################################################
### code chunk number 26: 18PlanesResidVFits
###################################################
Fits = fitted(Planes.aov2)
plot(Residuals~Fits)


###################################################
### code chunk number 27: 19PlanesLogResidQQNorm
###################################################
Residuals = residuals(Planes.aov3)
qqnorm(Residuals)
qqline(Residuals)


###################################################
### code chunk number 28: 19PlanesLogResidVFits
###################################################
Fits = fitted(Planes.aov3)
plot(Residuals~Fits)


###################################################
### code chunk number 29: Ch8.rnw:251-253
###################################################
summary(aov(Time~Gender*Lights, data=GoGoGo))
summary(aov(Time~Lights*Gender, data=GoGoGo))


###################################################
### code chunk number 30: Ch8.rnw:261-262
###################################################
bartlett.test(Price~Discipline, data=Textbooks)


