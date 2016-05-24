### R code from vignette source 'random-aids.Rnw'

###################################################
### code chunk number 1: random-aids.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: random-aids.Rnw:19-21 (eval = FALSE)
###################################################
## library(catdata)
## data(aids)


###################################################
### code chunk number 3: random-aids.Rnw:26-27 (eval = FALSE)
###################################################
## library(mgcv)


###################################################
### code chunk number 4: random-aids.Rnw:33-35 (eval = FALSE)
###################################################
## gammaids<-gamm(cd4 ~ s(time) + drugs + partners + s(cesd) + s(age), 
##                random=list(person=~1),  family=poisson(link=log), data=aids)


###################################################
### code chunk number 5: random-aids.Rnw:40-41 (eval = FALSE)
###################################################
## summary(gammaids$gam)


###################################################
### code chunk number 6: random-aids.Rnw:46-47 (eval = FALSE)
###################################################
## plot(gammaids$gam,ylab=" ",cex.lab=1.8,cex.axis=1.5,select=1)


###################################################
### code chunk number 7: random-aids.Rnw:50-51 (eval = FALSE)
###################################################
## plot(gammaids$gam,ylab=" ",cex.lab=1.8,cex.axis=1.5,select=2)


###################################################
### code chunk number 8: random-aids.Rnw:55-56 (eval = FALSE)
###################################################
## plot(gammaids$gam,ylab=" ",cex.lab=1.8,cex.axis=1.5,select=3)


###################################################
### code chunk number 9: random-aids.Rnw:58-59 (eval = FALSE)
###################################################
## detach(package:mgcv)


