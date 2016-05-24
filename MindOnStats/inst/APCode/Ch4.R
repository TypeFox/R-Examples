### R code from vignette source 'Ch4.rnw'

###################################################
### code chunk number 1: setup
###################################################
source("GenericSettings.R")


###################################################
### code chunk number 2: Ch4.rnw:15-16
###################################################
data(Bike, package="MindOnStats")


###################################################
### code chunk number 3: Ch4.rnw:19-21
###################################################
Bike2 = Bike[Bike$Type=="Bike", ]
str(Bike2)


###################################################
### code chunk number 4: 1Bike2StripPlotJitterWindows (eval = FALSE)
###################################################
## windows(7, 4)
## stripchart(Speed~Gender, data=Bike2, method="jitter", pch=20, xlab="Speed (km/h)")


###################################################
### code chunk number 5: 1Bike2StripPlotJitterNoWindowsByGender
###################################################
stripchart(Speed~Gender, data=Bike2, method="jitter", pch=20, xlab="Speed (km/h)")


###################################################
### code chunk number 6: 1Bike2StripPlotStackByGender
###################################################
stripchart(Speed~Gender, data=Bike2, method="stack", pch=20, xlab="Speed (km/h)")


###################################################
### code chunk number 7: 2Bike2BoxplotByGender
###################################################
boxplot(Speed~Gender, data=Bike2, ylab="Speed (km/h)", xlab="Gender", main="Speed of cyclists by gender")


###################################################
### code chunk number 8: Ch4.rnw:49-51
###################################################
data(Mobiles, package="MindOnStats")
str(Mobiles)


###################################################
### code chunk number 9: 3MobilesBoxplotByGenderPlan
###################################################
boxplot(Bill~paste(Gender,PlanType), data=Mobiles, ylab="Bill ($)", main="Monthly bill by gender and whether prepaid or plan", cex.axis=0.9)


###################################################
### code chunk number 10: Ch4.rnw:67-69
###################################################
data(Fish, package="MindOnStats")
str(Fish)


###################################################
### code chunk number 11: Ch4.rnw:74-76
###################################################
SortedLengths = sort(Fish$Length)
SortedLengths


###################################################
### code chunk number 12: Ch4.rnw:79-81
###################################################
mean(Fish$Length)
median(Fish$Length)


###################################################
### code chunk number 13: Ch4.rnw:88-90
###################################################
length(SortedLengths)
FishLengthsNoShark = SortedLengths[-58]


###################################################
### code chunk number 14: Ch4.rnw:95-99
###################################################
range(FishLengthsNoShark)
sd(FishLengthsNoShark)
IQR(FishLengthsNoShark)
var(FishLengthsNoShark)


###################################################
### code chunk number 15: Ch4.rnw:102-105
###################################################
fivenum(FishLengthsNoShark)
summary(FishLengthsNoShark)
quantile(FishLengthsNoShark, probs=c(0.25, 0.75))


###################################################
### code chunk number 16: Ch4.rnw:114-115
###################################################
lnBills = log(Mobiles$Bill)


###################################################
### code chunk number 17: Ch4.rnw:120-121
###################################################
quantile(FishLengthsNoShark, probs=0.90)


###################################################
### code chunk number 18: Ch4.rnw:124-126
###################################################
min(FishLengthsNoShark)
max(FishLengthsNoShark)


