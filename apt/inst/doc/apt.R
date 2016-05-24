### R code from vignette source 'apt.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: apt.Rnw:85-90
###################################################
library(apt)
data(daVich)
daVi <- y <- daVich[, 1] 
daCh <- x <- daVich[, 2]
bsStat(daVich)


###################################################
### code chunk number 2: apt.Rnw:95-98
###################################################
library(urca)
adf.xa <- ur.df(daCh,  type=c("trend"), lags=3)
adf.xa


###################################################
### code chunk number 3: apt.Rnw:102-103
###################################################
plot(adf.xa)


