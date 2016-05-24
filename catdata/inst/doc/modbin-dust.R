### R code from vignette source 'modbin-dust.Rnw'

###################################################
### code chunk number 1: modbin-dust.Rnw:11-14 (eval = FALSE)
###################################################
## library(catdata)
## data(dust)
## attach(dust)


###################################################
### code chunk number 2: modbin-dust.Rnw:17-19 (eval = FALSE)
###################################################
## dustlogitnon1=glm(bronch ~ dust+years, family=binomial, data=dust[(dust$smoke==0),])
## summary(dustlogitnon1)


###################################################
### code chunk number 3: modbin-dust.Rnw:22-25 (eval = FALSE)
###################################################
## dustlogitnon2 <- glm(bronch ~ dust+years, family=binomial, 
##                      data=dust[(dust$smoke==0)&(dust$dust<10),])
## summary(dustlogitnon2)


###################################################
### code chunk number 4: modbin-dust.Rnw:30-32 (eval = FALSE)
###################################################
## dustlogit1 <- glm(bronch ~ dust+years+smoke, family=binomial, data=dust)
## summary(dustlogit1)


###################################################
### code chunk number 5: modbin-dust.Rnw:35-38 (eval = FALSE)
###################################################
## dustlogit2 <- glm(bronch ~ dust+years+smoke, family=binomial, 
##                   data=dust[(dust$dust<20),])
## summary(dustlogit2)


