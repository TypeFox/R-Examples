### R code from vignette source 'pendensity.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: pendensity.Rnw:46-47
###################################################
library(pendensity)


###################################################
### code chunk number 2: pendensity.Rnw:52-55
###################################################
set.seed(27)
y <- rnorm(100)
test <- pendensity(y~1)


###################################################
### code chunk number 3: fig1
###################################################
plot(test)


###################################################
### code chunk number 4: pendensity.Rnw:76-80
###################################################
set.seed(27)
x <- rep(c(0,1),200)
y <- rnorm(400,x*0.5,1)
test2 <- pendensity(y~as.factor(x))


###################################################
### code chunk number 5: fig2
###################################################
plot(test2)


###################################################
### code chunk number 6: pendensity.Rnw:99-100
###################################################
test.equal(test2)


###################################################
### code chunk number 7: pendensity.Rnw:144-148
###################################################
x <- c(rep(0,50),rep(1,100),rep(2,100))
y <- rnorm(250,x,1)
x <- as.factor(x)
test3 <- pendensity(y~x)


###################################################
### code chunk number 8: fig11
###################################################
plot(test3,latt=TRUE)


###################################################
### code chunk number 9: pendensity.Rnw:155-156
###################################################
test.equal(test3)


###################################################
### code chunk number 10: pendensity.Rnw:172-174
###################################################
points <- c(1,-0.5,0,0.5,1)
plot(test,val=points)


###################################################
### code chunk number 11: pendensity.Rnw:179-181
###################################################
points <- c(1,-0.5,0,0.5,1)
plot(test2,val=points)


###################################################
### code chunk number 12: pendensity.Rnw:191-197
###################################################
data(Allianz)
form<-'%d.%m.%y'
time.Allianz <- strptime(Allianz[,1],form)
data.Allianz <- Allianz[which(time.Allianz$year==106),2]
d.Allianz <- diff(data.Allianz)
density.Allianz <- pendensity(d.Allianz~1)


###################################################
### code chunk number 13: fig3
###################################################
plot(density.Allianz)


###################################################
### code chunk number 14: fig4
###################################################
plot(density.Allianz,latt=TRUE)


###################################################
### code chunk number 15: pendensity.Rnw:222-225
###################################################
data.Allianz <- Allianz[which(time.Allianz$year==106|time.Allianz$year==107),2]
d.Allianz <- diff(data.Allianz)
density.Allianz2 <- pendensity(d.Allianz~as.factor(time.Allianz$year))


###################################################
### code chunk number 16: fig5
###################################################
plot(density.Allianz2)


###################################################
### code chunk number 17: fig6
###################################################
plot(density.Allianz2,latt=TRUE)


###################################################
### code chunk number 18: pendensity.Rnw:246-247
###################################################
test.equal(density.Allianz2)


###################################################
### code chunk number 19: pendensity.Rnw:254-260
###################################################
set.seed(27)
y <- rnorm(400)
density1 <- pendensity(y~1,no.base=10)
density2 <- pendensity(y~1,no.base=15)
density3 <- pendensity(y~1,no.base=20)
density4 <- pendensity(y~1,no.base=25)              


###################################################
### code chunk number 20: fig7
###################################################
plot(density1)


###################################################
### code chunk number 21: fig8
###################################################
plot(density2)


###################################################
### code chunk number 22: fig9
###################################################
plot(density3)


###################################################
### code chunk number 23: fig10
###################################################
plot(density4)


