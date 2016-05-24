### R code from vignette source 'binary-vaso.Rnw'

###################################################
### code chunk number 1: binary-vaso.Rnw:11-14 (eval = FALSE)
###################################################
## library(catdata)
## data(vaso)
## attach(vaso)


###################################################
### code chunk number 2: binary-vaso.Rnw:18-20 (eval = FALSE)
###################################################
## y <- vaso$vaso
## y[vaso$vaso==2] <- 0


###################################################
### code chunk number 3: binary-vaso.Rnw:23-25 (eval = FALSE)
###################################################
## vaso1 <- glm(y ~ vol + rate, family=binomial)
## summary(vaso1)


###################################################
### code chunk number 4: binary-vaso.Rnw:28-30 (eval = FALSE)
###################################################
## vaso2 <- glm(y ~ I(exp(vol)) + I(exp(rate)), family=binomial)
## summary(vaso2)


