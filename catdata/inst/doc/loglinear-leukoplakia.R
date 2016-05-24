### R code from vignette source 'loglinear-leukoplakia.Rnw'

###################################################
### code chunk number 1: loglinear-leukoplakia.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=60)


###################################################
### code chunk number 2: loglinear-leukoplakia.Rnw:19-23 (eval = FALSE)
###################################################
## library(catdata)
## data(leukoplakia)
## library(vcdExtra)
## leukoplakia<-expand.dft(leukoplakia,freq="Freq")


###################################################
### code chunk number 3: loglinear-leukoplakia.Rnw:28-30 (eval = FALSE)
###################################################
## mosaicplot(~ Leukoplakia+Smoker+Alcohol,main ="", data = leukoplakia, 
##            color = TRUE, cex.axis=1.0)


###################################################
### code chunk number 4: loglinear-leukoplakia.Rnw:33-35 (eval = FALSE)
###################################################
## mosaicplot(~ Leukoplakia+Smoker,main ="", data = leukoplakia, color = TRUE, 
##            cex.axis=1.0)


###################################################
### code chunk number 5: loglinear-leukoplakia.Rnw:38-40 (eval = FALSE)
###################################################
## detach(package:vcdExtra)
## detach(package:gnm)


