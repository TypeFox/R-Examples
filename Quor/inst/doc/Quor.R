### R code from vignette source 'Quor.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: quorinst (eval = FALSE)
###################################################
## install.packages("Quor_VERSION.tar.gz",
##                  repos=NULL,type="source")


###################################################
### code chunk number 2: <quor
###################################################
library("Quor")


###################################################
### code chunk number 3: example2
###################################################
data(gleason7)
d <- list(x1 = gleason7[1:5,1], x2 = gleason7[,2])
conf.statement(d,verbose=FALSE)


