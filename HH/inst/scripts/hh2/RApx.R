### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/RApx.tex'

###################################################
### code chunk number 1: RApx.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: RApx.tex:49-66
###################################################
hhcode("install.packagesLM.R", '
install.packages(c("HH","RcmdrPlugin.HH","RcmdrPlugin.mosaic",
                   "fortunes","ggplot2","shiny","gridExtra",
                   "gridBase","Rmpfr","png","XLConnect",
                   "matrixcalc", "sem", "relimp", "lmtest",
                   "markdown", "knitr", "effects", "aplpack",
                   "RODBC", "TeachingDemos",
                   "gridGraphics", "gridSVG"),
                 dependencies=TRUE)

## This is the sufficient list (as of 16 August 2015) of packages
## needed in order to install the HH package.  Should
## additional dependencies be declared by any of these packages
## after that date, the first use of "library(HH)" after the
## installation might ask for permission to install some more
## packages.
')


###################################################
### code chunk number 3: RApx.tex:97-102
###################################################
hhcode("install.packages.R", '
## Tell Windows that R should have the same access to the
## outside internet that is granted to Internet Explorer.
## setInternet2()  ## this line is defunct beginning with R_3.3.0
')


###################################################
### code chunk number 4: RApx.tex:104-108
###################################################
hhcode("install.RcmdrW.R", '
install.packages("Rcmdr",
                 dependencies=TRUE)
')


###################################################
### code chunk number 5: RApx.tex:271-276
###################################################
## hhcapture("simpleR.Rout", '
## Simple R session
3 + 4
pnorm(c(-1.96, -1.645, -0.6745, 0, 0.6745, 1.645, 1.96))
## ')


###################################################
### code chunk number 6: RApx.tex:304-307
###################################################
## hhcode("manual.R", '
system.file("../../doc/manual")
## ')


###################################################
### code chunk number 7: RApx.tex:309-312
###################################################
## hhcode("manualWindows.R", '
WindowsPath(system.file("../../doc/manual"))
## ')


###################################################
### code chunk number 8: RApx.tex:361-364
###################################################
## hhcode("HHscript2.R", '
HHscriptnames()
## ')


###################################################
### code chunk number 9: RApx.tex:366-369
###################################################
## hhcode("HHscript1.R", '
HHscriptnames(edition=1)
## ')


###################################################
### code chunk number 10: RApx.tex:371-374
###################################################
## hhcode("HHscriptW2.R", '
WindowsPath(HHscriptnames())
## ')


###################################################
### code chunk number 11: RApx.tex:376-379
###################################################
## hhcode("HHscriptW1.R", '
WindowsPath(HHscriptnames(edition=1))
## ')


