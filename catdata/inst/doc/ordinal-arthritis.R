### R code from vignette source 'ordinal-arthritis.Rnw'

###################################################
### code chunk number 1: ordinal-arthritis.Rnw:12-13 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: ordinal-arthritis.Rnw:19-23 (eval = FALSE)
###################################################
## arthritis <- data.frame(drug=c(rep("new agent", 24+37+21+19+6), 
## rep("active control", 11+51+22+21+7)), 
## assessment=c(rep(1,24), rep(2,37), rep(3,21), rep(4,19), rep(5,6), rep(1,11), 
## rep(2,51), rep(3,22), rep(4,21), rep(5,7)))


###################################################
### code chunk number 3: ordinal-arthritis.Rnw:26-27 (eval = FALSE)
###################################################
## library(VGAM)


###################################################
### code chunk number 4: ordinal-arthritis.Rnw:32-34 (eval = FALSE)
###################################################
## cumart <- vglm(assessment ~ drug, family=cumulative(parallel=TRUE, link="logit"), 
##               data=arthritis)


###################################################
### code chunk number 5: ordinal-arthritis.Rnw:37-38 (eval = FALSE)
###################################################
## summary(cumart)


###################################################
### code chunk number 6: ordinal-arthritis.Rnw:41-42 (eval = FALSE)
###################################################
## detach(package:VGAM)


