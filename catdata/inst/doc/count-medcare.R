### R code from vignette source 'count-medcare.Rnw'

###################################################
### code chunk number 1: count-medcare.Rnw:11-14 (eval = FALSE)
###################################################
## library(catdata)
## data(medcare)
## attach(medcare)


###################################################
### code chunk number 2: count-medcare.Rnw:18-21 (eval = FALSE)
###################################################
## med1=glm(ofp ~ hosp+healthpoor+healthexcellent+numchron+age+married+school,
##          family=poisson,data=medcare[male==1 & ofp<=30,])
## summary(med1)


###################################################
### code chunk number 3: count-medcare.Rnw:24-27 (eval = FALSE)
###################################################
## med2=glm(ofp ~ hosp+healthpoor+healthexcellent+numchron+age+married+school,
##          family=quasipoisson,data=medcare[male==1 & ofp<=30,])
## summary(med2)


###################################################
### code chunk number 4: count-medcare.Rnw:31-35 (eval = FALSE)
###################################################
## library(MASS)
## med3=glm.nb(ofp ~ hosp+healthpoor+healthexcellent+numchron+age+married+school,
##             data=medcare[male==1 & ofp<=30,])
## summary(med3)


###################################################
### code chunk number 5: count-medcare.Rnw:38-38 (eval = FALSE)
###################################################
## 


###################################################
### code chunk number 6: count-medcare.Rnw:41-42 (eval = FALSE)
###################################################
## library(pscl)


###################################################
### code chunk number 7: count-medcare.Rnw:44-47 (eval = FALSE)
###################################################
## med4=zeroinfl(ofp ~ hosp+healthpoor+healthexcellent+numchron+age+married+school|1,
##               data=medcare[male==1 & ofp<=30,])
## summary(med4)


###################################################
### code chunk number 8: count-medcare.Rnw:50-53 (eval = FALSE)
###################################################
## med5=zeroinfl(ofp ~ hosp+healthpoor+healthexcellent+numchron+age+married+school,
##               data=medcare[male==1 & ofp<=30,])
## summary(med5)


###################################################
### code chunk number 9: count-medcare.Rnw:56-59 (eval = FALSE)
###################################################
## med6=hurdle(ofp ~ hosp+healthpoor+healthexcellent+numchron+age+married+school|1
##             ,data=medcare[male==1 & ofp<=30,])
## summary(med6)


###################################################
### code chunk number 10: count-medcare.Rnw:61-64 (eval = FALSE)
###################################################
## med7=hurdle(ofp ~ hosp+healthpoor+healthexcellent+numchron+age+married+school,
##             data=medcare[male==1 & ofp<=30,])
## summary(med7)


