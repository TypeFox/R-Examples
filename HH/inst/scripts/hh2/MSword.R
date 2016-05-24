### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/MSword.tex'

###################################################
### code chunk number 1: MSword.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: MSword.tex:150-154
###################################################
## hhcapture("RVar.Rout", '
k <- c(0, 1, 2, 15, 16, 17, 26, 27, 28)
cbind(k=k, var=apply(cbind(10^k + 1, 10^k + 2, 10^k + 3), 1, var))
## ')


###################################################
### code chunk number 3: MSword.tex:175-179
###################################################
## hhcapture("ExcelError.Rout", '
sprintf("%+13.13a", 28334198897217900000000)
as.matrix(sprintf("%+13.13a", c(10^27 + 1, 10^27 + 2, 10^27 + 3)))
## ')


