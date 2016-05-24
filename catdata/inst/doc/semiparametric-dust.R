### R code from vignette source 'semiparametric-dust.Rnw'

###################################################
### code chunk number 1: semiparametric-dust.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: semiparametric-dust.Rnw:19-22 (eval = FALSE)
###################################################
## library(mgcv)
## library(catdata)
## data(dust)


###################################################
### code chunk number 3: semiparametric-dust.Rnw:28-32 (eval = FALSE)
###################################################
## gamdust1 <- gam(bronch ~ s(dust,years), family=binomial, 
##                 data=dust[(dust$dust<10) & (dust$smoke==1),])
## 
## plot(gamdust1, pers=TRUE)


###################################################
### code chunk number 4: semiparametric-dust.Rnw:37-39 (eval = FALSE)
###################################################
## gamdust2<- gam(bronch ~ s(dust) + s(years), family=binomial, 
##                data=dust[(dust$dust<10) & (dust$smoke==1),])


###################################################
### code chunk number 5: semiparametric-dust.Rnw:41-42 (eval = FALSE)
###################################################
## plot(gamdust2, select=1)


###################################################
### code chunk number 6: semiparametric-dust.Rnw:45-46 (eval = FALSE)
###################################################
## plot(gamdust2, select=2)


