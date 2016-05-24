### R code from vignette source 'multivariate-birth1.Rnw'

###################################################
### code chunk number 1: multivariate-birth1.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: multivariate-birth1.Rnw:20-23 (eval = FALSE)
###################################################
## library(catdata)
## data(birth)
## attach(birth)


###################################################
### code chunk number 3: multivariate-birth1.Rnw:29-32 (eval = FALSE)
###################################################
## intensive <- rep(0,length(Intensive))
## intensive[Intensive>0] <- 1
## Intensive <- intensive


###################################################
### code chunk number 4: multivariate-birth1.Rnw:37-40 (eval = FALSE)
###################################################
## previous <- Previous
## previous[previous>1] <- 2 
## Previous <- previous


###################################################
### code chunk number 5: multivariate-birth1.Rnw:43-44 (eval = FALSE)
###################################################
## library(VGAM)


###################################################
### code chunk number 6: multivariate-birth1.Rnw:49-52 (eval = FALSE)
###################################################
## Birth <- as.data.frame(na.omit(cbind(Intensive, Cesarean, Sex, Weight, Previous, 
## AgeMother)))
## detach(birth)


###################################################
### code chunk number 7: multivariate-birth1.Rnw:57-61 (eval = FALSE)
###################################################
## bivarlogit <- vglm(cbind(Intensive , Cesarean) ~ as.factor(Sex) + Weight + 
## as.factor(Previous) + AgeMother, binom2.or(zero=NULL), data=Birth)
## 
## summary(bivarlogit)


###################################################
### code chunk number 8: multivariate-birth1.Rnw:64-65 (eval = FALSE)
###################################################
## detach(package:VGAM)


