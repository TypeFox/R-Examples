### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/RExcelApx.tex'

###################################################
### code chunk number 1: RExcelApx.tex:8-9
###################################################
library(HH)


###################################################
### code chunk number 2: RExcelApx.tex:66-77
###################################################
hhcode("install.packages.RExcel.R", '
## Tell Windows that R should have the same access to the outside
## internet that is granted to Internet Explorer.
## setInternet2()  ## this line is defunct beginning with R_3.3.0

install.packages(c("rscproxy","rcom"),
                 repos="http://rcom.univie.ac.at/download",
                 type="binary",
                 lib=.Library)
library(rcom)
comRegisterRegistry()
')


###################################################
### code chunk number 3: RExcelApx.tex:259-282
###################################################
if (FALSE) {
## The data in this figure is included with the installer
##    RthroughExcelWorkbooksInstaller_1.2-10.exe
## and is therefore available only on Windows through the RExcel
## Add-In button.  It is not available on other operating systems.
##
## We did some sleight of hand on this figure.  We changed the colors to
## a more intense blue and red.  We changed the order of the levels of the
## Gender factor to make the order in the key match the order in the graph.
## We reversed the colors.

StudentData$Gender <-
   factor(StudentData$Gender, levels=c("male","female"))

trellis.par.set(list(superpose.symbol=list(col=
  likertColor(10, colorFunctionOption = "default")[c(10,2)]
)))

xyplot(Shoesize ~ Size + SizeFather + SizeMother, layout=c(1, 3),
  groups=Gender, type="p", pch=16, auto.key=list(border=TRUE),
  par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'),
  y=list(relation='same')), data=StudentData)
}


