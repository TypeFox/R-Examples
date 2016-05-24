### R code from vignette source 'Ch6.rnw'

###################################################
### code chunk number 1: setup
###################################################
source("GenericSettings.R")


###################################################
### code chunk number 2: Ch6.rnw:31-34
###################################################
x=1:5
x
cumsum(x)


###################################################
### code chunk number 3: Ch6.rnw:41-42
###################################################
dbinom(7,  size=10, prob=0.487)


###################################################
### code chunk number 4: Ch6.rnw:47-50
###################################################
x=10:15
dbinom(x, size=15, prob=0.5)
sum(dbinom(x, size=15, prob=0.5))


###################################################
### code chunk number 5: Ch6.rnw:53-54
###################################################
pbinom(9, size=15, prob=0.5, lower.tail=FALSE) 


###################################################
### code chunk number 6: Ch6.rnw:59-60
###################################################
pnorm(25, mean=26.42, sd=4.3)


###################################################
### code chunk number 7: Ch6.rnw:63-64
###################################################
pnorm(30, mean=26.42, sd=4.3, lower.tail=FALSE)


###################################################
### code chunk number 8: Ch6.rnw:67-69
###################################################
pnorm(c(25,30), mean=26.42, sd=4.3)
diff(pnorm(c(25,30), mean=26.42, sd=4.3))


###################################################
### code chunk number 9: Ch6.rnw:80-82
###################################################
data(Bike, package="MindOnStats")
MaleCyclists=Bike$Speed[Bike$Gender=="Male"&Bike$Type=="Bike"]


###################################################
### code chunk number 10: 12MaleCyclistsQQNorm
###################################################
qqnorm(MaleCyclists)
qqline(MaleCyclists)


###################################################
### code chunk number 11: Ch6.rnw:95-97
###################################################
data(TimePerception, package="MindOnStats")
str(TimePerception)


###################################################
### code chunk number 12: 13FiveSecHist
###################################################
attach(TimePerception)
hist(FiveSec)


###################################################
### code chunk number 13: 13FiveSecQQNorm
###################################################
qqnorm(FiveSec)
qqline(FiveSec)
detach(TimePerception)


###################################################
### code chunk number 14: 14TenSecHist
###################################################
attach(TimePerception)
hist(TenSec)


###################################################
### code chunk number 15: 14TenSecQQNorm
###################################################
qqnorm(TenSec)
qqline(TenSec)
detach(TimePerception)


