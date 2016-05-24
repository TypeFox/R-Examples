### R code from vignette source 'count-encephalitis.Rnw'

###################################################
### code chunk number 1: count-encephalitis.Rnw:11-14 (eval = FALSE)
###################################################
## library(catdata)
## data(encephalitis)
## attach(encephalitis)


###################################################
### code chunk number 2: count-encephalitis.Rnw:19-22 (eval = FALSE)
###################################################
## BAV <- country
## BAV[BAV==2] <-0
## TIME <- year


###################################################
### code chunk number 3: count-encephalitis.Rnw:27-29 (eval = FALSE)
###################################################
## enc1 <- glm(count ~ TIME+I(TIME^2)+BAV+TIME*BAV, family = poisson)
## summary(enc1)


###################################################
### code chunk number 4: count-encephalitis.Rnw:33-35 (eval = FALSE)
###################################################
## enc2 <- glm(count ~ TIME+I(TIME^2)+BAV+TIME*BAV, family = gaussian("identity"))
## summary(enc2)


###################################################
### code chunk number 5: count-encephalitis.Rnw:39-42 (eval = FALSE)
###################################################
## enc3 <- glm(count ~ TIME+I(TIME^2)+BAV+TIME*BAV, family = gaussian("log"), 
##             start=enc1$coef)
## summary(enc3)


